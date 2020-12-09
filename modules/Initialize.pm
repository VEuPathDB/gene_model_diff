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

=cut

=head1 CONTACT
	
	Please email comments or questions to help@veupathdb.org
	
=cut

=head1 NAME

Initialize

=head1 SYNOPSIS

	use Initialize;
	
	Initialize::load_gene_set($dbh,$config,$validation_file,$source,$gff_file,$fasta_file,$dsn,$user,$pass);

=head1 DESCRIPTION

This module loads the GFF into a SeqFeature Store database either in memory or disk if connection details is given.
The genes that passes the validation is inserted into the gene model table.

=head1 Author

	Mikkel B Christensen

=head1 METHODS

=cut

package Initialize;
use strict;
use warnings;
use DBI;
use Data::Dumper;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use Digest::MD5 qw(md5_hex);
use Validate;

use Bio::Seq;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;



my $verbose = 0;
my $db;
my $loader;

=head2 load_gene_set

 Title: load_gene_set
 Usage: Initialize::load_gene_set($dbh,$config,$validation_file,$source,$gff_file,$fasta_file,$dsn,$user,$pass)
 Function: load the gff into a SeqFeature Store database, then loads validated genes into the gene_model table.
 Returns: Counts of genes: before load, Unfinished, faild validation, total loaded.
 Args: Database handle object, config object,source of GFF, GFF file,FASTA file,database connection details (dns),user,pass
=cut

sub load_gene_set {
	my($dbh,$config,$validation_file,$source,$gff_file,$fasta_file,$dsn,$user,$pass) = @_;
	my $not_finished=0;
	my $not_validated=0;
	my $total_loaded=0;
	my $infoLog = Log::Log4perl->get_logger("infologger");
	my $pre_loaded = _gff_load($gff_file,$fasta_file,$dsn,$user,$pass);
	
	my @genes = $db->get_features_by_type('gene');
    open my $validation_fh,'>>',$validation_file;
	foreach my $gene (@genes){
		
		my %attb = $gene->attributes;
	
		if($source eq 'cap'){
			unless(defined($attb{status}->[0]) and ($attb{status}->[0] eq 'Finished')){
				$not_finished++;
				next;
			}			
		}
		
		my $passed_validation = Validate::validate_gene($gene,$config,$validation_fh);
		
		if($source eq 'cap' and $passed_validation < 0){
			$not_validated++;
			$infoLog->info("Gene $attb{load_id}->[0] did not paas validation");
		}
				
		my ($gene_model,$CDS_present) = _build_gene_model($gene);
	    
		if(%{$gene_model}){
			_insert_gene_model($gene_model,$dbh,$source);
			_insert_exon($gene_model,$dbh,$source);
			if($CDS_present){
				_insert_CDS($gene_model,$dbh,$source);
			}
		}
		
		$total_loaded++;
	}
	
	return($pre_loaded,$not_finished,$not_validated,$total_loaded);
}

sub _gff_load {
	my($gff_file_name,$fasta_file_name,$dsn,$user,$pass) = @_;
	warn "load GFF";
	$db=();
	$loader=();
	if($dsn){
		$db = Bio::DB::SeqFeature::Store->new(
			-adaptor	=> 'DBI::mysql',
             -dsn	=> $dsn,
             -user   => $user,
             -pass   => $pass,
             -create	=> 1 
		);
	}else{
		$db = Bio::DB::SeqFeature::Store->new(
			-adaptor => 'memory',
		);
	}
  
  	$loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(
  		-store   => $db,
  		-verbose => $verbose,
  	);
  	
  	open my $gff_fh,'<', $gff_file_name;
  	my $loaded;
  	if($fasta_file_name){
  		open my $fasta_fh,'<',$fasta_file_name;
  		$loaded = $loader->load($gff_fh,$fasta_fh);  		
  	}else{my $loaded = $loader->load($gff_fh);}
  	
  	return $loaded;  	
}

sub _build_gene_model {
	my($gene) = @_;	
	my %attb = $gene->attributes;
	my %gene_model;
	my $CDS_present = 0;
	$gene_model{gene_id} = $attb{load_id}->[0];
	$gene_model{validation_error_code} = $attb{validation_error_code}->[0];
	
	my @mRNAs = $gene->get_SeqFeatures('mRNA');
	foreach my $mRNA (@mRNAs){
		my %mRNA_model;
		my %rna_attb = $mRNA->attributes;
		my $mRNA_validation_error_code = $rna_attb{validation_error_code}->[0];
		
		$mRNA_model{mRNA_id} = $rna_attb{load_id}->[0];
		$mRNA_model{validation_error_code} = $mRNA_validation_error_code;
		
		
		my @CDS = $mRNA->get_SeqFeatures('CDS');
		if(scalar(@CDS) > 0){
			$CDS_present = 1;
			my $CDS_sequence = get_CDS($mRNA,'');
			
			my $prot_seq     = get_translation($CDS_sequence,$mRNA->strand);
			my $md5_checksum = md5_hex($prot_seq);
			
			my %cds_attb  = $CDS[0]->attributes;
			$mRNA_model{cds_start}        = $CDS[0]->start;
			$mRNA_model{cds_end}          = $CDS[0]->end;
			$mRNA_model{CDS_Parent_id}    = $cds_attb{parent_id}->[0];
			$mRNA_model{CDS_md5_checksum} = $md5_checksum;
			$mRNA_model{cds_error_code}   = $mRNA_validation_error_code;
		}
		
		
		my @exons = $mRNA->get_SeqFeatures('exon');
				
		foreach my $exon (@exons){
		    my %exon_attb = $exon->attributes;
			my %exon_model;
			
			$exon_model{exon_id}   = $exon_attb{load_id}->[0]?$exon_attb{load_id}->[0]:$exon_attb{parent_id}->[0] . $exon_attb{rank}->[0];#make uniq exon ID
			$exon_model{scaffold}  = $exon->seq_id;
			$exon_model{strand}    = $exon->strand;
			$exon_model{start}     = $exon->start;
			$exon_model{end}       = $exon->end;
			
			push @{$mRNA_model{exon}}, \%exon_model;
		}
		push@{$gene_model{mRNA}},\%mRNA_model;	
	}
	
	return (\%gene_model,$CDS_present);
}



sub get_CDS {
	my($mrna_obj,$validation_fh) = @_;
	my $infoLog = Log::Log4perl->get_logger("infologger");
	my $seq_id = $mrna_obj->seq_id;
	my $strand = $mrna_obj->strand;
	
	my @CDS = $mrna_obj->get_SeqFeatures('CDS');
	my %attb = $CDS[0]->attributes;
	#only expects one CDS per mRNA, code will break if more.
	if(scalar @CDS > 1){croak("more than one CDS for mRNA,
							   only expects one CDS per mRNA, code will break if more. $attb{load_id}->[0]");}
	my $CDS_start = $CDS[0]->start;
	my $CDS_end   = $CDS[0]->end;
	my $CDS_sequence;
	
	my @exons  = $mrna_obj->get_SeqFeatures('exon');
	@exons = sort {$a->start <=> $b->start} @exons;
	
	#need to splice the mRNA
	if(scalar @exons == 1){
		$CDS_sequence = get_sequence($seq_id,$CDS_start,$CDS_end);
	}elsif(scalar @exons > 1){
		my $exon_start = shift @exons;
		my $exon_end   = pop @exons;
		my $exon_start_endpos = $exon_start->end;
		my $exon_end_startpos = $exon_end->start;
		my $cds_in_one_exon = 0;	
		
		#do not include 5' UTR
		while($exon_start_endpos < $CDS_start){		
			if(!scalar @exons){
				$cds_in_one_exon=1;
				last;
			}
			$exon_start = shift @exons;
			$exon_start_endpos = $exon_start->end;
		}
		
		#do not include 3' UTR
		while( $CDS_end < $exon_end_startpos){			
			if(!scalar @exons){
				$cds_in_one_exon=1;
				last
			}
			$exon_end   = pop @exons;						
			$exon_end_startpos = $exon_end->start;						
		}
		
		if($cds_in_one_exon){
			$CDS_sequence = get_sequence($seq_id,$CDS_start,$CDS_end);
		}else{				
			my $CDS_start_seq = get_sequence($seq_id,$CDS_start,$exon_start_endpos);
			
			my $CDS_end_seq = get_sequence($seq_id,$exon_end_startpos,$CDS_end); 
			
			chomp($CDS_start_seq,$CDS_end_seq);
		   
			if(scalar(@exons)){
				my $exon_sequence;
				foreach my $exon (@exons){	    		
					my $exon_start = $exon->start;
					my $exon_end   = $exon->end;
					$exon_sequence .= get_sequence($seq_id,$exon_start,$exon_end);	    		
					#print $validation_fh "middle exon\t" . $exon_sequence . "\n"; 
				}
				chomp($CDS_start_seq,$exon_sequence,$CDS_end_seq);
				$CDS_sequence = $CDS_start_seq . $exon_sequence . $CDS_end_seq;
			}else{$CDS_sequence = $CDS_start_seq . $CDS_end_seq;}		
		}
	}
	unless($CDS_sequence){
		$infoLog->info("NO CDS sequence was returned for mRNA, $attb{load_id}->[0]");	
	}
	return $CDS_sequence;	
}

sub get_translation {
	my($CDS_sequence,$strand) = @_;
	
	my $seqobj = Bio::Seq->new(-seq => $CDS_sequence);
	if($strand eq '-' or $strand eq '-1'){
			$seqobj = $seqobj->revcom();
	}
	my $prot_seq = $seqobj->translate->seq;
	
	return $prot_seq;
}



sub get_sequence {
	my($scaffold,$start,$end) = @_;
	my $seq = $db->fetch_sequence($scaffold,$start,$end);
	return $seq;	
}

sub _insert_exon{ 
	my ($hash,$dbh,$source) = @_;
	
	my $insert_sql = "insert exon(exon_id,
								  source,
								  scaffold,
								  strand,
								  start,
								  end
					 )
					 select ?,?,?,?,?,?;";
	my $sth	= $dbh->prepare($insert_sql);						
	foreach my $mRNA (@{$hash->{mRNA}}){
		foreach my $exon (@{$mRNA->{exon}}){
				$sth->bind_param(1,$exon->{exon_id});
				$sth->bind_param(2,$source);
				$sth->bind_param(3,$exon->{scaffold});
				$sth->bind_param(4,$exon->{strand});
				$sth->bind_param(5,$exon->{start});
				$sth->bind_param(6,$exon->{end});
				
				$sth->execute();
		}
	}
	
		
}

sub _insert_CDS{
	my ($hash,$dbh,$source) = @_;
	
	my $insert_sql ="insert cds(cds_parent_id,
								start,
								end,
								md5_checksum,
								cds_error_code
					 )
					 select ?,?,?,?,?;";
	
	my $sth	= $dbh->prepare($insert_sql);
	foreach my $mRNA (@{$hash->{mRNA}}){
		
		$sth->bind_param(1,$mRNA->{CDS_Parent_id});
		$sth->bind_param(2,$mRNA->{cds_start});
		$sth->bind_param(3,$mRNA->{cds_end});
		$sth->bind_param(4,$mRNA->{CDS_md5_checksum});
		$sth->bind_param(5,$mRNA->{cds_error_code});
				
		$sth->execute();	
	}
}


sub _insert_gene_model{
	my ($hash,$dbh,$source) = @_;
	my $insert_sql = "insert gene_model(
										exon_id,
										transcript_id,
										gene_id,
										source,
										error_code
									   )
							select ?,?,?,?,?;";
	my $sth	= $dbh->prepare($insert_sql);
	my $gene_id = $hash->{gene_id};
	foreach my $mRNA (@{$hash->{mRNA}}){
		foreach my $exon (@{$mRNA->{exon}}){
				
				$sth->bind_param(1,$exon->{exon_id});
				$sth->bind_param(2,$mRNA->{mRNA_id});
				$sth->bind_param(3,$hash->{gene_id});
				$sth->bind_param(4,$source);
				$sth->bind_param(5,$hash->{validation_error_code});
				
				$sth->execute();
		}
	}	
}

1;
