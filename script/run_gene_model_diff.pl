=head1 LICENSE

Copyright [2017-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 NAME

gene_model_diff.pl - Merging two genesets

=head1 SYNOPSIS

gene_model_diff.pl [options]

	Options:
		--config		server and directory specification file
		--speciesFile		list of species text file

=head1 OPTIONS

=over 4

=item B<--config>

Configuration file with database connectivity details and working directories specifications

=item B<--speciesFile>

Text file with list of species name that are also the name of the subdirectories with the GFF and Fasta files are located

=back

=head1 DESCRIPTION

B<This program> will read a fasta file and two GFF files (the original and and updated geneset) and produces a GFF file with the merged geneset and a file with genes ids to delete from the original geneset.

=cut

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use DBI;
use Config::IniFiles;
use Getopt::Long;
use Pod::Usage;

use FindBin;
use lib("$FindBin::Bin/../modules/.");
use ExonMapping;
use GeneClusters;
use TranscriptLinks;
use TranscriptMapping;
use GeneMapping;
use AnnotationEvents;
use Initialize;


my $event_new_count;
my $event_change_count;
my $event_merge_count;
my $event_split_count;
my $delete_total_count;
my $load_total_count;
my $gff_gene_count;

my $identical_gene_count;
my $new_gene_count;
my $changed_gene_count;
my @merged_gene_counts;
my @split_gene_counts;

my $dbh;

#----------------Get options---------------------#
my %options;
my $result = GetOptions(\%options,
		'config|f=s',
		'speciesFile=s') || pod2usage(2);
#------------------------------------------------#

my $config = readConfig($options{config});
my $datadir= $config->val('Data','datadir');
my $scriptdir = $FindBin::Bin;
my $species_list = get_species($options{speciesFile});
my $log_file = $config->val('Data','log_file');
Log::Log4perl->init($log_file);

foreach my $species (@{$species_list}){
	chomp($species);
	run_species($species,1);
}

sub run_species {
	my($species,$run_count) = @_;
	
	$event_new_count =0;
	$event_change_count=0;
	$event_merge_count=0;
	$event_split_count=0;
	$delete_total_count=0;
	$load_total_count=0;
	$gff_gene_count=0;

	$identical_gene_count=0;
	$new_gene_count=0;
	$changed_gene_count=0;
	@merged_gene_counts=0;
	@split_gene_counts=0;
	
	my $cap_total_loaded=0;
	
	warn "Running $species\n";
	
    create_database($config,$species);		
	$cap_total_loaded = load_data_base($config,$species);
	
	$dbh = get_dbh($config,$species);	
	if($cap_total_loaded){
		run_gene_set_comparison($dbh);
		
		
		get_events($dbh);
		write_events($config,$species);
		write_delete_list($config,$species);
		write_gff_to_load($config,$species);
		write_summary_counts($config,$species);
	}else{warn "No genes was loaded for $species\n"}	
}

sub load_data_base {
	my($config,$species) = @_;
	chomp($config,$species);
	warn "Running load database for $species\n";
	
		
	my $datadir= $config->val('Data','datadir');	
	my $host = $config->val('Database','host');
	my $port = $config->val('Database','port');
	my $user = $config->val('Database','user');
	my $pass = $config->val('Database','pass');
	my $load_database = $config->val('Database','load_database');
	
	my $cap_gff    = "$datadir/$species/cap.gff";
	my $core_gff   = "$datadir/$species/core.gff";
	my $cap_fasta  = "$datadir/$species/cap.fasta";
	my $core_fasta = "$datadir/$species/core.fasta";
	my $validation_file = "$datadir/$species/gff_validation_error.txt";
	
	my $dns  = "dbi:mysql:$load_database:$host:$port";	
	$dbh = get_dbh($config,$species);
	
	my $pruned_core_gff = prune_gff_by_scaffold($cap_gff,$core_gff);
	
	my($cap_pre_loaded,$cap_not_finished,$cap_not_validated,$cap_total_loaded) = Initialize::load_gene_set($dbh,$config,$validation_file,'cap',$cap_gff,$cap_fasta,$dns,$user,$pass);
	warn "Stats: loaded=$cap_pre_loaded, not_finished=$cap_not_finished, not_validated=$cap_not_validated, total_loaded=$cap_total_loaded\n";
	if($cap_total_loaded){
		my($vb_pre_loaded,$vb_not_finished,$vb_not_validated,$vb_total_loaded) = Initialize::load_gene_set($dbh,$config,$validation_file,'vb',$pruned_core_gff,$core_fasta,$dns,$user,$pass);
	}
	return $cap_total_loaded;
 	
}

sub prune_gff_by_scaffold {
	my($cap_gff,$core_gff) = @_;
	my %cap_scaffold;
	
	my $pruned_core_gff = $core_gff . '.pruned';
	open my $cap_fh,'<',$cap_gff;
	open my $core_fh,'<',$core_gff;
	open my $pruned_core_gff_fh,'>',$pruned_core_gff;
	
	while (my $line = <$cap_fh>) {
		chomp $line;
		next if $line =~ /^#/;
		my $scaffold = (split/\t/,$line)[0];
		next if not $scaffold;
		$cap_scaffold{$scaffold} = 1;
	}
	
	while (my $line = <$core_fh>) {		
		next if $line =~ /^###/;
		if ($line =~ /^#/) {
			print $pruned_core_gff_fh $line;
		} else {
			my $scaffold = (split/\t/,$line)[0];
      next if not $scaffold;
			if (exists $cap_scaffold{$scaffold}) {
				print $pruned_core_gff_fh $line;
			}
		}
	}
	return $pruned_core_gff;	
}

sub run_gene_set_comparison {
	my($dbh) = @_;

	ExonMapping::work_out_exon_mapping($dbh);
	GeneClusters::work_out_gene_clusters($dbh);
	GeneClusters::calculate_cluster_summary($dbh);
	TranscriptLinks::work_out_transcript_links($dbh);	
	TranscriptMapping::resolve_transcript_mappings($dbh);
	GeneMapping::resolve_maptype_cluster($dbh);

}

sub get_events {
	my($dbh) = @_;
	$identical_gene_count = AnnotationEvents::get_identical_gene($dbh);
 	$new_gene_count =       AnnotationEvents::get_new_gene($dbh);
 	$changed_gene_count =   AnnotationEvents::get_changed_genes($dbh);
 	@merged_gene_counts =   AnnotationEvents::get_merge($dbh);
 	@split_gene_counts =    AnnotationEvents::get_splits($dbh);
}

sub write_events {
	my($config,$species) = @_;
	chomp($config,$species);
	my %duplicate_ids;
	my $datadir= $config->val('Data','datadir');
	my $gff_file_dir = "$datadir/$species";
	open my $file_handle,'>',"$gff_file_dir/annotation_events.txt";
	my $sql = 'select vb_gene_id,cap_gene_id,events from gene_events;';
	
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	
	my $events_ref = $sth->fetchall_arrayref;
	if(scalar @{$events_ref} < 1){
		#die "No events found\n;";
	}
	
	foreach my $row (@{$events_ref}){
		my($vb_gene_id,$cap_gene_id,$event) = @{$row};
		if((defined($vb_gene_id) and exists $duplicate_ids{$vb_gene_id}) or (defined($cap_gene_id) and exists $duplicate_ids{$cap_gene_id})){
			print 	$file_handle "Duplicate\t$vb_gene_id~$cap_gene_id\n//\n";
		}elsif($event eq 'new_gene'){
			$duplicate_ids{$cap_gene_id} =1;
		}else{
			$duplicate_ids{$vb_gene_id} =1;
			$duplicate_ids{$cap_gene_id} =1
		}
		
		if($event eq 'identical'){
			print $file_handle "Ge\t$vb_gene_id~$cap_gene_id\n//\n";
		}
		
		
		if($event eq 'change_gene'){
			print $file_handle "Ge\t$vb_gene_id=$cap_gene_id\n//\n";
			$event_change_count++;
		}elsif($event eq 'new_gene'){
			print $file_handle "Ge\t+$cap_gene_id\n//\n";
			$event_new_count++;
		}elsif($event eq "merge_gene"){
			print $file_handle "Ge\t$vb_gene_id>$cap_gene_id\n//\n";
			$event_merge_count++;
		}elsif($event eq "split_gene"){
			print $file_handle "Ge\t$vb_gene_id<$cap_gene_id\n//\n";
			$event_split_count++;
		}
		
				
	}
	
}

sub write_delete_list {
	my($config,$species) = @_;
	chomp($config,$species);
	my $datadir= $config->val('Data','datadir');
	
	my $sql = "select distinct(vb_gene_id) from gene_mappings where map_type <> 'new' and map_type <> 'identical';";
	my $array_ref = $dbh->selectall_arrayref($sql);
	
	my $gff_file_dir = "$datadir/$species";
	open my $file_handle,'>',"$gff_file_dir/old_genes_delete_from_core.txt";
	
	foreach my $id (@{$array_ref}){				
		print $file_handle "$id->[0]\n";
		$delete_total_count++;	
	}
	
}

sub write_gff_to_load {
	my($config,$species) = @_;
	chomp($config,$species);
	my %loaded_hash;
	my %parents_hash;
	my $datadir= $config->val('Data','datadir');
	my $gff_file_dir= "$datadir/$species";
	
	my $sql = "select distinct(cap_gene_id) from gene_mappings where map_type <> 'identical';";
	my $array_ref = $dbh->selectall_arrayref($sql);
	
	open my $file_handle,'>',"$gff_file_dir/loaded_genes_delete_from_WA.txt";
	
	for my $id (@{$array_ref}){
		print $file_handle "$id->[0]\n";
		$loaded_hash{$id->[0]}=1;
		$load_total_count++;
	}
	
	my $cap_gff= "$gff_file_dir/cap.gff";
	open my $cap_gff_fh,'<', $cap_gff or die "can't open $cap_gff\n";
	open my $load_gff_fh,'>',"$gff_file_dir/genes_2_load.gff" or die "can't open gene_2_load.gff\n";
	print $load_gff_fh "##gff-version 3\n";
	
	while( my $line = <$cap_gff_fh>){
		chomp $line;
		
		my @colunms = split/\t/,$line;
		next unless($colunms[2]);
		if($colunms[2] eq 'gene'){
			my($ID) = $colunms[8]=~/ID=([\w\-]+)\;/;
			if(exists $loaded_hash{$ID}){
				print $load_gff_fh "$colunms[0]\tVectorBase\t$colunms[2]\t$colunms[3]\t$colunms[4]\t$colunms[5]\t$colunms[6]\t$colunms[7]\t$colunms[8]\n";
				$gff_gene_count++;
			}
		}elsif($colunms[2] eq 'mRNA'){			
			my($parent) = $colunms[8]=~/Parent=([\w\-]+)\;/;
			my($ID) = $colunms[8]=~/ID=([\w\-]+)\;/;
			if(exists $loaded_hash{$parent}){
				$parents_hash{$ID} = 1;
				print $load_gff_fh "$colunms[0]\tVectorBase\t$colunms[2]\t$colunms[3]\t$colunms[4]\t$colunms[5]\t$colunms[6]\t$colunms[7]\t$colunms[8]\n";
			}		
		}elsif($colunms[2] eq 'exon' or $colunms[2] eq 'CDS'){
			my($parent) = $colunms[8]=~/Parent=([\w\-]+)\;/;
			if(exists $parents_hash{$parent}){
				print $load_gff_fh "$colunms[0]\tVectorBase\t$colunms[2]\t$colunms[3]\t$colunms[4]\t$colunms[5]\t$colunms[6]\t$colunms[7]\t$colunms[8]\n";
			}					
		}	
	}
}

sub write_summary_counts {
	my($config,$species) = @_;
    chomp($config,$species);
	my $datadir= $config->val('Data','datadir');
	my $gff_file_dir= "$datadir/$species";
	
   open  my $file_fh,'>',"$gff_file_dir/summary_counts.txt";	
   print $file_fh  "Identical genes $identical_gene_count\n";
   print $file_fh  "New gene events:\t$event_new_count ($new_gene_count)\n";
   print $file_fh  "Changed gene events:\t$event_change_count ($changed_gene_count)\n";
   print $file_fh  "Merge gene events:\t$event_merge_count ($merged_gene_counts[0],$merged_gene_counts[1],$merged_gene_counts[2])\n";
   print $file_fh  "Split gene events:\t$event_split_count ($split_gene_counts[0],$split_gene_counts[1],$split_gene_counts[2])\n";
   print $file_fh  "Total genes deleted:\t$delete_total_count\n";
   print $file_fh  "Total genes loaded:\t$load_total_count\n";
   print $file_fh  "GFF gene count:\t$gff_gene_count\n";		
}

sub get_species {
	my($speciesFile) = @_;
	
	unless($speciesFile){ croak("No species file supplied");}
	unless(-e $speciesFile){ croak("species file $speciesFile do not exists");}
	
	open my $file_h,'<',$speciesFile;
	
	my @species = <$file_h>;
	
	return \@species;
}

sub readConfig{
    my ($configFile) = @_;
    unless($configFile){ croak("No config file supplied");}
    unless(-e $configFile){ croak("config file $configFile do not exists");}
    
    $config = Config::IniFiles->new(-file => $configFile);

    return $config;
}

sub create_database {
	my($config,$species) = @_;
	chomp($config,$species); 
	my $host = $config->val('Database','host');
	my $port = $config->val('Database','port');
	my $user = $config->val('Database','user');
	my $pass = $config->val('Database','pass');
	my $prefix= $config->val('Database','prefix');
	my $database_name = $config->val('Database','database');
		
	my $database = $prefix .'_' . $species . '_'  . $database_name;
	
	my $dns  = "dbi:mysql:;host=$host;port=$port";
	my $dbh=DBI->connect( $dns,$user, $pass,{RaiseError => 1}) or die "cant connect as $user";
	
	$dbh->do("drop database if exists $database;");
	$dbh->do("create database $database;");
	
	`mysql --host=$host --port=$port --user=$user  --password=$pass $database  < $scriptdir/create_classification_database.sql`;
	croak("Failed to load table definitions in create_classification_database.sql to database $database on host $host") if $?;
}

sub get_dbh {
	my($config,$species) = @_;
	chomp($config,$species); 
	my $host = $config->val('Database','host');
	my $port = $config->val('Database','port');
	my $user = $config->val('Database','user');
	my $pass = $config->val('Database','pass');
	my $prefix= $config->val('Database','prefix');
	my $database_name = $config->val('Database','database');
		
	my $database = $prefix .'_' . $species . '_'  . $database_name;
	
	my $dns  = "dbi:mysql:$database:$host:$port";
	my $dbh=DBI->connect( $dns,$user, $pass,{RaiseError => 1})
		or die "cant connect to $database as $user";
	
	return $dbh;
}




