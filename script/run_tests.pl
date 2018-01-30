use strict;
use warnings;
use Test::More qw/no_plan/;
use Carp;
use Data::Dumper;
use DBI;
use lib('../modules/.');


use ExonMapping;
use TranscriptLinks;
use TranscriptMapping;
use GeneMapping;
use AnnotationEvents;
use Initialize;
use GeneModel;
use GeneClusters;
use Config::IniFiles;
use Getopt::Long;

#----------------Get options---------------------#
my %options;
my $result = GetOptions(\%options,
		'config|f=s',
		'testFile=s');
#------------------------------------------------#

my $config       = readConfig($options{config});
my $datadir      = $config->val('Database','datadir');
my $test_list = get_species($options{testFile});#test list
my $dbh;
my $testFile;
my $species_list;

my %actions=(
	no_change=>\&no_change,
	no_change_isoforms=>\&no_change_isoforms,
	utr_change=>\&utr_change,
	new_gene=>\&new_gene,
	added_exon=>\&added_exon,
	remove_exon=>\&remove_exon,
	alter_exon=>\&alter_exon,
	gene_in_intron=>\&gene_in_intron,
	iso_form=>\&iso_form,
	iso_form_lost=>\&iso_form_lost,
	gene_merge=>\&gene_merge,
	gene_split=>\&gene_split,
	multi_gene_events=>\&multi_gene_events,
	gene_cluster=>\&gene_cluster,
	gene_cluster_merge_split=>\&gene_cluster_merge_split,
	gene_validation_ok=>\&gene_validation_ok,
	gene_validation_3exons_ok=>\&gene_validation_3exons_ok,
	gene_validation_no_start=>\&gene_validation_no_start,
	gene_validation_no_stop=>\&gene_validation_no_stop,
	interface_test_1=>\&interface_test_1,
	interface_test_2=>\&interface_test_2,
	interface_test_3=>\&interface_test_3,
	three_tests_dbs=>\&three_tests_dbs
);

`perl run_gene_model_diff.pl --config $options{config}  --species $options{testFile}`;
foreach my $test (@{$test_list}){
		chomp $test;
		print "------RUNNING $test --------\n"; 
		$actions{$test}->($test);		
}



sub get_events {
	my $sql = 'select vb_gene_id,cap_gene_id,events from gene_events;';
	my $sth = $dbh->prepare($sql);
	$sth->execute();
	
	my $events_ref = $sth->fetchall_arrayref;
	if(scalar @{$events_ref} < 1){
		my @array = ['faild','faild'];
		$events_ref = \@array;	
	}
	
	return $events_ref;
	
}
#######Test##############
sub utr_change{
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	my $events_ref = get_events();
		
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'GB0001'and $cap_gene_id eq 'G0001' and $events eq 'change_gene'), 'UTR change was indetified correct');
	}
}

sub gene_validation_no_stop{
	my($test) = @_;
	$dbh = get_dbh($config,$test);
	my $results_ref = GeneModel::get_distinct_id_by_source($dbh,'gene','cap');	
	print Dumper($results_ref);
	ok(scalar @{$results_ref} == 0, 'Gene whitout stop codon not loaded');
	
}

sub gene_validation_no_start{
	my($test) = @_;
	$dbh = get_dbh($config,$test);	
	my $results_ref = GeneModel::get_distinct_id_by_source($dbh,'gene','cap');	
	print Dumper($results_ref);
	ok(scalar @{$results_ref} == 0, 'Gene whitout start codon not loaded');
	
}

sub gene_validation_ok{	
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	my $events_ref = get_events();
		
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq '08135d5a-14fe-4f94-a37e-7ee7b4947bf3'and $cap_gene_id eq '08135d5a-14fe-4f94-a37e-7ee7b4947bf3' and $events eq 'identical'), 'validation of genes were correct');
	}	
}

sub gene_validation_3exons_ok{
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	my $events_ref = get_events();
		
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq '2809cfc5-804b-4ac2-abff-d9bc3d230b2b'and $cap_gene_id eq '2809cfc5-804b-4ac2-abff-d9bc3d230b2b' and $events eq 'identical'), 'validation of genes were correct');
	}		
}

sub no_change {
    my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	my $events_ref = get_events();
		
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'GB0001'and $cap_gene_id eq 'G0001' and $events eq 'identical'), 'No change identical genes were indetified correct');
	}	
}

sub no_change_isoforms {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	my $events_ref = get_events();
		
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0001' and $events eq 'identical'), 'No change isoform identical genes were indetified correct');
	}
	
}

sub new_gene {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	my $events_ref = get_events();
	
	my($vb_gene_id,$cap_gene_id,$events);
	my @test_values = ('','G0002','new_gene','new gene indetified correct','GB0001','G0001','identical','identical genes were indetified correct in new gene');
	my $test_count =0;
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		next if($events eq 'identical');
		print $vb_gene_id . $test_values[$test_count] . $cap_gene_id . $test_values[$test_count + 1] . $events . $test_values[$test_count + 2] . "\n";
		ok(($vb_gene_id eq $test_values[$test_count] and $cap_gene_id eq $test_values[$test_count + 1] and  $events eq $test_values[$test_count + 2]),$test_values[$test_count + 3]);	
		$test_count = 4;		
	}
	
}

sub added_exon {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);

    my $events_ref = get_events();
	
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0001' and $events eq 'change_gene'), 'Added exon change_gene genes were indetified correct');
	}
	
}

sub remove_exon {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
    my $events_ref = get_events();
	
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0001' and $events eq 'change_gene'), 'Lost exon change_gene genes were indetified correct');
	}	
}

sub alter_exon {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
    my $events_ref = get_events();
	
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0001' and $events eq 'change_gene'), 'Alter Exon change_gene genes were indetified correct');
	}		
}

sub gene_in_intron {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	my $events_ref = get_events();
	
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		next if($events eq 'identical');
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq ''and $cap_gene_id eq 'G0002' and $events eq 'new_gene'), 'new_gene gene in intron were indetified correct');
	}		
}

sub iso_form {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	my $events_ref = get_events();
	
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0001' and $events eq 'change_gene'), 'new_isoform were indetified correct');
	}		
}

sub iso_form_lost {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	my $events_ref = get_events();
	
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0001' and $events eq 'change_gene'), 'lost_isoform were indetified correct');
	}		
}



sub gene_merge {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);	
	
	my $events_ref = get_events();
	
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok((($vb_gene_id eq 'G0002:G0001' or $vb_gene_id eq 'G0001:G0002') and $cap_gene_id eq 'G0001' and $events eq 'merge_gene'), 'Gene merge were indetified correct');
	}			
}

sub gene_split {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	my $events_ref = get_events();
	
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0002:G0001' and $events eq 'split_gene'), 'split gene were indetified correct');
	}			
}

sub multi_gene_events {
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	
	my $events_ref = get_events();
	foreach my $row (@{$events_ref}){
		my($vb_gene_id,$cap_gene_id,$events) = @{$row};	
		if($cap_gene_id eq 'G0001'){
			ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0001' and $events eq 'identical'), 'No change identical genes were indetified correct');
		}elsif($cap_gene_id eq 'G0002'){
			ok(($vb_gene_id eq 'G0002'and $cap_gene_id eq 'G0002' and $events eq 'change_gene'), 'Alter Exon change_gene genes were indetified correct');		
		}elsif($vb_gene_id eq 'G0003'){
			ok(($vb_gene_id eq 'G0003'and $cap_gene_id eq 'G0004:G0003' and $events eq 'split_gene'), 'split gene were indetified correct');
		}elsif($cap_gene_id eq 'G0005'){
			ok(($vb_gene_id eq '' and $cap_gene_id eq 'G0005' and  $events eq 'new_gene'),'new gene indetified correct');	
		}elsif($cap_gene_id eq 'G0006'){
			ok((($vb_gene_id eq 'G0007:G0006' or $vb_gene_id eq 'G0006:G0007') and $cap_gene_id eq 'G0006' and $events eq 'merge_gene'), 'Gene merge were indetified correct');
		}
	}		
}

sub gene_cluster{
	my($test) = @_;	
	$dbh = get_dbh($config,$test);
	
	GeneClusters::work_out_gene_clusters($dbh);
	my $clusters_ref = GeneClusters::get_gene_clusters($dbh);
	ok(scalar @{$clusters_ref} == 3 , 'There is 3 genes in clusters');
	foreach my $row (@{$clusters_ref}){
		my($cluster_id,$gene_id,$source) = @{$row};
		print "$cluster_id,$gene_id,$source\n";
		if($cluster_id == 1 and $source eq 'cap'){ok($gene_id eq 'G0001', 'correct cap gene');}
		if($cluster_id == 2 and $source eq 'cap'){ok($gene_id eq 'G0002', 'correct cap gene');}
		if($cluster_id == 1 and $source eq 'vb'){ok($gene_id eq 'GB0001','correct vb gene');}			
	}
	
	GeneClusters::calculate_cluster_summary($dbh);
	my $cluster_summary = GeneClusters::get_cluster_summary($dbh);
	my($gene_cluster_id,$cap_gene_count,$cap_trans_count,$vb_gene_count,$vb_trans_count) = @{$cluster_summary->[0]};
	if($gene_cluster_id == 1){ok($cap_gene_count == 1 and $cap_trans_count ==1 and $vb_gene_count == 1 and $vb_trans_count == 1,"cluster summary ok");}
			
}

sub gene_cluster_merge_split {
	my($test) = @_;
	$dbh = get_dbh($config,$test);
	GeneClusters::work_out_gene_clusters($dbh);
	my $clusters_ref = GeneClusters::get_gene_clusters($dbh);
	ok(scalar @{$clusters_ref} == 6 , 'There is 6 genes in cluster');
	
		
}

sub readConfig{
    my ($configFile) = @_;
    unless($configFile){ die "No config file supplied\n";}
    unless(-e $configFile){ die "config file $configFile do not exists\n";}
    
    $config = Config::IniFiles->new(-file => $configFile);

    return $config;
}

sub get_species {
	my($speciesFile) = @_;
	
	unless($speciesFile){ die "No species file supplied\n";}
	unless(-e $speciesFile){ die "species file $speciesFile do not exists\n";}
	
	open my $file_h,'<',$speciesFile;
	
	my @species = <$file_h>;
	
	return \@species;
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
	my $dbh=DBI->connect( $dns,"ensrw", $pass,{RaiseError => 1})
		or die "can't connect to $database as $user";
	
	return $dbh;
}

sub get_dbh_old {
	my($database,$database_server,$pass) = @_;
	
	my $dbh='';
	if($database_server eq 'mysql'){
		my $database= 'mbc_' . $database;
		my $host    ='localhost';
		my $dns  = "dbi:mysql:$database:$host";
		   $dbh=DBI->connect( $dns,"mikkel", $pass,{RaiseError => 1})
			   or die "can't connect to $database as ensrw";
		
	}elsif($database_server eq 'mysql-devel-1'){
		my $database= 'mbc_' . $database;
		my $host    ='mysql-eg-devel-1.ebi.ac.uk';
		my $dns  = "dbi:mysql:$database:$host:4126";
		   $dbh=DBI->connect( $dns,"ensrw", $pass,{RaiseError => 1})
			or die "can't connect to $database as ensrw";
	}else{
		my $database=$database;
		my $host    ='mysql-eg-prod-vb.ebi.ac.uk';
		my $dns  = "dbi:mysql:$database:$host:4440";
		   $dbh=DBI->connect( $dns,"ensrw", $pass,{RaiseError => 1})
			   or die "can't connect to $database as ensrw";
	}

	return $dbh;
}

sub three_tests_dbs{
	`perl ../run_classify_annotation_events.pl --config $options{config}  --species $testFile`;	
	
	$dbh = get_dbh($config,$species_list->[0]);#no_change
	my $events_ref = get_events();
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'GB0001'and $cap_gene_id eq 'G0001' and $events eq 'identical'), 'No change identical genes were indetified correct');
	}
	
	$dbh = get_dbh($config,$species_list->[1]);#no_change_isoforms
	$events_ref = get_events();
	$vb_gene_id=$cap_gene_id=$events='';
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0001' and $events eq 'identical'), 'No change isoform identical genes were indetified correct');
	}
	
	$dbh = get_dbh($config,$species_list->[2]);#new_gene
	$events_ref = get_events();
	$vb_gene_id=$cap_gene_id=$events='';
	my @test_values = ('GB0001','G0001','identical','identical genes were indetified correct in new gene','','G0002','new_gene','new gene indetified correct');
	my $test_count =0;
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		
		print $vb_gene_id . $test_values[$test_count] . $cap_gene_id . $test_values[$test_count + 1] . $events . $test_values[$test_count + 2] . "\n";
		ok(($vb_gene_id eq $test_values[$test_count] and $cap_gene_id eq $test_values[$test_count + 1] and  $events eq $test_values[$test_count + 2]),$test_values[$test_count + 3]);	
		$test_count = 4;		
	}
	
}

sub interface_test_1{
	`perl run_gene_model_diff.pl --config $options{config}  --species $options{speciesFile}`;
		
	$dbh = get_dbh($config,$species_list->[0]);
	my $events_ref = get_events();
		
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq 'GB0001'and $cap_gene_id eq 'G0001' and $events eq 'identical'), 'No change identical genes were indetified correct');
	}
	
}

sub interface_test_2{
	`perl run_gene_model_diff.pl --config $options{config}  --species $options{speciesFile}`;
	
	$dbh = get_dbh($config,$species_list->[0]);
	my $events_ref = get_events();
	foreach my $row (@{$events_ref}){
		my($vb_gene_id,$cap_gene_id,$events) = @{$row};	
		if($cap_gene_id eq 'G0001'){
			ok(($vb_gene_id eq 'G0001'and $cap_gene_id eq 'G0001' and $events eq 'identical'), 'No change identical genes were indetified correct');
		}elsif($cap_gene_id eq 'G0002'){
			ok(($vb_gene_id eq 'G0002'and $cap_gene_id eq 'G0002' and $events eq 'change_gene'), 'Alter Exon change_gene genes were indetified correct');		
		}elsif($vb_gene_id eq 'G0003'){
			ok(($vb_gene_id eq 'G0003'and $cap_gene_id eq 'G0004:G0003' and $events eq 'split_gene'), 'split gene were indetified correct');
		}elsif($cap_gene_id eq 'G0005'){
			ok(($vb_gene_id eq '' and $cap_gene_id eq 'G0005' and  $events eq 'new_gene'),'new gene indetified correct');	
		}elsif($cap_gene_id eq 'G0006'){
			ok((($vb_gene_id eq 'G0007:G0006' or $vb_gene_id eq 'G0006:G0007') and $cap_gene_id eq 'G0006' and $events eq 'merge_gene'), 'Gene merge were indetified correct');
		}
	}
}

sub interface_test_3{
	
	`perl run_gene_model_diff.pl --config $options{config}  --species $options{speciesFile}`;
	$dbh = get_dbh($config,$species_list->[0]);
	my $events_ref = get_events();
		
	my($vb_gene_id,$cap_gene_id,$events); 
	foreach my $row (@{$events_ref}){
		
		$vb_gene_id  = $row->[0];
		$cap_gene_id = $row->[1];
		$events      = $row->[2];
		print "$vb_gene_id:$cap_gene_id:$events\n";
		ok(($vb_gene_id eq '08135d5a-14fe-4f94-a37e-7ee7b4947bf3'and $cap_gene_id eq '08135d5a-14fe-4f94-a37e-7ee7b4947bf3' and $events eq 'identical'), 'validation of genes were correct');
	}

}