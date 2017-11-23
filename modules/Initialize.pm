=head1 LICENSE

Copyright [2017] EMBL-European Bioinformatics Institute

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

package Initialize;
use strict;
use warnings;
use DBI;
use Data::Dumper;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;

use Bio::Seq;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;



my $verbose = 0;
my $do_validation = 1;
my $db;
my $loader;

sub load_gene_set {
	my($dbh,$validation_file,$source,$gff_file,$fasta_file) = @_;
	my $not_finished=0;
	my $not_validated=0;
	my $total_loaded=0;
	my $pre_loaded = gff_load($gff_file,$fasta_file);
	
	my @genes = $db->get_features_by_type('gene');
    open my $validation_fh,'>',$validation_file if $source eq 'cap';
	foreach my $gene (@genes){
		my %attb = $gene->attributes;
		if($do_validation and $source eq 'cap'){
			#warn "gene: $attb{status}->[0]\n";
			unless(defined($attb{status}->[0]) and ($attb{status}->[0] eq 'Finished annotating')){
				$not_finished++;
				next;
			}			
		}
		
		if($do_validation and $source eq 'cap' and !validate_gene($gene,$validation_fh)){
			#warn "gene:$attb{load_id}->[0] not loaded\n";
			$not_validated++;
			next;
		}
		
		my %gene_model = build_gene_model($gene);
	    
		if(%gene_model){
			insert_exon(\%gene_model,$dbh,$source);
			insert_CDS(\%gene_model,$dbh,$source);
			insert_gene_model(\%gene_model,$dbh,$source);	
		}
		
		$total_loaded++;
	}
	
	return($pre_loaded,$not_finished,$not_validated,$total_loaded);
}

sub build_gene_model {
	my($gene) = @_;	
	my %attb = $gene->attributes;
	my %gene_model;
	$gene_model{gene_id} = $attb{load_id}->[0];
	
	my @mRNAs = $gene->get_SeqFeatures('mRNA');
	foreach my $mRNA (@mRNAs){
		my %mRNA_model;
		my %rna_attb = $mRNA->attributes;
		$mRNA_model{mRNA_id} = $rna_attb{load_id}->[0];
		
		my @CDS = $mRNA->get_SeqFeatures('CDS');
		
		my %cds_attb  = $CDS[0]->attributes;
		$mRNA_model{cds_start}     = $CDS[0]->start;
		$mRNA_model{cds_end}       = $CDS[0]->end;
		$mRNA_model{CDS_Parent_id} = $cds_attb{parent_id}->[0];
				
		my @exons = $mRNA->get_SeqFeatures('exon');
				
		foreach my $exon (@exons){
		    my %exon_attb = $exon->attributes;
			my %exon_model;
			#make uniq ex
			$exon_model{exon_id}   = $exon_attb{load_id}->[0]?$exon_attb{load_id}->[0]:$exon_attb{parent_id}->[0] . $exon_attb{rank}->[0];#make uniq exon ID
			$exon_model{scaffold}  = $exon->seq_id;
			$exon_model{strand}    = $exon->strand;
			$exon_model{start}     = $exon->start;
			$exon_model{end}       = $exon->end;
			
			push @{$mRNA_model{exon}}, \%exon_model;
		}
		push@{$gene_model{mRNA}},\%mRNA_model;	
	}
	
	return %gene_model;
}

sub validate_gene {
	my($gene,$validation_fh) = @_;
	my %gene_attb = $gene->attributes;
	
	
	my %gene_hash;
	$gene_hash{scaffold} = $gene->seq_id;
	$gene_hash{strand}   = $gene->strand;
	$gene_hash{start}    = $gene->start;
	$gene_hash{end}      = $gene->end;
	
	my @RNAs = $gene->get_SeqFeatures('mRNA');
	my $passed_validation = 1;
	my $validation_string;
	
	foreach my $mRNA (@RNAs){
		my %mRNA_hash;
		my %mRNA_attb = $mRNA->attributes;
		if((exists $mRNA_attb{'no-ATG'}) or (exists $mRNA_attb{'no-STOP'})){
						$validation_string = "partial prediction";
						$passed_validation = 0;
						last;
		}
		
		eval{				
				my $mRNA_ID   = $mRNA_attb{load_id}->[0];
				$mRNA_hash{scaffold} = $mRNA->seq_id;
				$mRNA_hash{strand}   = $mRNA->strand;
				$mRNA_hash{start}    = $mRNA->start;
				$mRNA_hash{end}      = $mRNA->end;
				
				my $validation_text;
				if(!check_subfeature(\$validation_text,$mRNA_ID,\%gene_hash,\%mRNA_hash,'Gene:mRNA')){
					$validation_string .= "mRNA:$mRNA_ID|$validation_text;";
					$passed_validation = 0;
				}
				
				my @exons = $mRNA->get_SeqFeatures('exon');
				
				foreach my $exon (@exons){
					my %exon_hash;
					my %exon_attb = $exon->attributes;
					my $exon_ID = $exon_attb{load_id}->[0];
					$exon_hash{scaffold} = $exon->seq_id;
					$exon_hash{strand}   = $exon->strand;
					$exon_hash{start}    = $exon->start;
					$exon_hash{end}      = $exon->end;
				
					if(!check_subfeature(\$validation_text,$exon_ID,\%mRNA_hash,\%exon_hash,'mRNA:exon')){
						$validation_string .=  "exon:$exon_ID|$validation_text;";
						$passed_validation = 0;
					}
				
				}
				
				my $CDS_sequence = get_CDS($mRNA,$validation_fh);
				my $prot_seq     = get_translation($CDS_sequence,$mRNA_hash{strand});
				
				if(!check_seq($prot_seq,\$validation_text)){										
					$validation_string .= "mRNA:$mRNA_ID";
					$validation_string .= "$validation_text;";
					$passed_validation = 0;
				}
		};
		if($@){
			print $validation_fh $@;
			$passed_validation = 0;
		}
	}
	
	if(!$passed_validation){
		my $gene_id = $gene_attb{load_id}->[0];
		my $owner =   $gene_attb{owner}->[0];
		
		print $validation_fh "gene:$gene_id\t$owner\t$validation_string\n";
	
	}
	
	return $passed_validation;
}

sub check_subfeature{
	my($validation_text,$ID,$feature,$sub_feature,$key) = @_;
	my $passed = 1;	
	my ($feature_name,$sub_feature_name) = split /:/,$key;
	$$validation_text='';
	if($feature->{scaffold} ne $sub_feature->{scaffold}){
		$$validation_text .= "error:Scaffold not indentical|proof:[$feature_name scaffold $feature->{scaffold}, $sub_feature_name scaffold $sub_feature->{scaffold}];";
		$passed = 0;
	}
	
	if($feature->{strand} ne $sub_feature->{strand}){
		$$validation_text .= "error:Strand not indentical|proof:[$feature_name strand $feature->{strand}, $sub_feature_name strand $sub_feature->{strand}];";
		$passed = 0;
	}
	
	if($feature->{start} > $sub_feature->{start}){
		$$validation_text .= "error:sub feature start is not within borders|proof:[$feature_name start $feature->{start}, $sub_feature_name start $sub_feature->{start}];";
		$passed = 0;
	}
	
	if($feature->{end} < $sub_feature->{end}){
		$$validation_text .= "error:sub feature end is not within borders|proof:[$feature_name end $feature->{end}, $sub_feature_name end $sub_feature->{end}];";
		$passed = 0;
	}
	return $passed;
}

sub gff_load {
	my($gff_file_name,$fasta_file_name) = @_;
	warn "load GFF";
	$db=();
	$loader=();
	$db = Bio::DB::SeqFeature::Store->new(
		-adaptor => 'memory',
    );
  
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

sub get_CDS {
	my($mrna_obj,$validation_fh) = @_;
	
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
		
		while($exon_start_endpos < $CDS_start){		
			if(!scalar @exons){
				$cds_in_one_exon=1;
				last;
			}
			$exon_start = shift @exons;
			$exon_start_endpos = $exon_start->end;
		}
			
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

sub check_seq {
	my($prot_seq,$validation_text) = @_;
	my $start_aa = substr($prot_seq,0,1);
	my $end_aa   = substr($prot_seq,-1,1);
	my $start_flag=0;
	my $stop_flag =0;
	if($start_aa eq 'M'){$start_flag++;};
	if($end_aa eq '*'){$stop_flag++;};
	my $internal_stop_count = -1;
	while($prot_seq=~m/\*/g){$internal_stop_count++;}
	$$validation_text = '';
	if(!$start_flag){
		$$validation_text .= "|error:No start codon";
	}
	if(!$stop_flag){
		$$validation_text .= "|error:No stop codon";
	}
	if($internal_stop_count > 0){
		$$validation_text .= "|error:Internal stop codon";
	}
	
	if($$validation_text){
		$$validation_text .= "|proof:$prot_seq;";	
	}
	
	if($start_flag and $stop_flag and ($internal_stop_count == 0)){
		return 1;
	}else{return 0;}
}

sub get_sequence {
	my($scaffold,$start,$end) = @_;
	my $seq = $db->fetch_sequence($scaffold,$start,$end);
	return $seq;	
}

sub insert_exon{ 
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

sub insert_CDS{
	my ($hash,$dbh,$source) = @_;
	
	my $insert_sql ="insert cds(cds_parent_id,
								start,
								end
					 )
					 select ?,?,?;";
	
	my $sth	= $dbh->prepare($insert_sql);
	foreach my $mRNA (@{$hash->{mRNA}}){
		
		$sth->bind_param(1,$mRNA->{CDS_Parent_id});
		$sth->bind_param(2,$mRNA->{cds_start});
		$sth->bind_param(3,$mRNA->{cds_end});
				
		$sth->execute();	
	}
}


sub insert_gene_model{
	my ($hash,$dbh,$source) = @_;
	my $insert_sql = "insert gene_model(
										exon_id,
										transcript_id,
										gene_id,
										source
									   )
							select ?,?,?,?;";
	my $sth	= $dbh->prepare($insert_sql);
	my $gene_id = $hash->{gene_id};
	foreach my $mRNA (@{$hash->{mRNA}}){
		foreach my $exon (@{$mRNA->{exon}}){
				
				$sth->bind_param(1,$exon->{exon_id});
				$sth->bind_param(2,$mRNA->{mRNA_id});
				$sth->bind_param(3,$hash->{gene_id});
				$sth->bind_param(4,$source);
				
				$sth->execute();
		}
	}	
}

sub set_do_validation {
	my($boolean)= @_;
	$do_validation = $boolean;
}

1;
