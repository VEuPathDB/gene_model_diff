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

Validate

=head1 SYNOPSIS

	use Validate;
	
	Validate::validate_gene($gene,$config,$validation_fh);

=head1 DESCRIPTION

This module validates the gene model. The following checks are applied: All sub-feature
must be on the same strand and within its parents boundaries. Coding sequence must have Start and Stop codon
or have a 'no-ATG'/'no-stop' tag to pass.

=head1 Author

	Mikkel B Christensen

=head1 METHODS

=cut

package Validate;

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use Initialize;

use Bio::Seq;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;


=head2 validate_gene

 Title: validate_gene
 Usage: Validate::validate_gene($gene,$config,$validation_fh)
 Function: validates the gene model
 Returns: validation/error code, which are given in the config file. 
 Args: SeqFeature gene object, config object, validation file handle
=cut

sub validate_gene {
	my($gene,$config,$validation_fh) = @_;
	my %gene_attb = $gene->attributes;
	my $errorLog = Log::Log4perl->get_logger("errorlogger");
	
	
	my %validation_error_code = ( 'GFF_mRNA'  => $config->val('Validate','GFF_mRNA'),
						   		 'GFF_exon'  => $config->val('Validate','GFF_exon'),
						   		 'CDS_start' => $config->val('Validate','CDS_start'),
						   		 'CDS_stop'  => $config->val('Validate','CDS_stop'),
						   		 'CDS_internal_stop' => $config->val('Validate','CDS_internal_stop'),
						   		 'runtime_error' => $config->val('Validate','runtime_error')
						   	   );
	my %validation_approved_email = ( 'approved_email' => $config->val('Validate','approved_email'));
	
	my %gene_hash;
	$gene_hash{scaffold} = $gene->seq_id;
	$gene_hash{strand}   = $gene->strand;
	$gene_hash{start}    = $gene->start;
	$gene_hash{end}      = $gene->end;
	
	my @RNAs = $gene->get_SeqFeatures('mRNA');
	my $gene_validation_status = 0;
	my $validation_string;
	
	foreach my $mRNA (@RNAs){
		my %mRNA_hash;
		my %mRNA_attb = $mRNA->attributes;
		my $mRNA_validation_status = 0;
		
		my $NO_ATG = 0;
		my $NO_STOP = 0;

		if(exists $mRNA_attb{'no-ATG'}){
			if($mRNA_attb{owner}->[0] eq $validation_approved_email{approved_email}){
				$NO_ATG = 2;
			}else{ $NO_ATG = 1;}
		}
		
		if(exists $mRNA_attb{'no-STOP'}){
			if($mRNA_attb{owner}->[0] eq $validation_approved_email{approved_email}){
				$NO_STOP = 2;
			}else{ $NO_STOP = 1;}
		}
				
		eval{				
				my $mRNA_ID   = $mRNA_attb{load_id}->[0];
				$mRNA_hash{scaffold} = $mRNA->seq_id;
				$mRNA_hash{strand}   = $mRNA->strand;
				$mRNA_hash{start}    = $mRNA->start;
				$mRNA_hash{end}      = $mRNA->end;
				
				my $validation_text;
				if(!_check_subfeature(\$validation_text,$mRNA_ID,\%gene_hash,\%mRNA_hash,'Gene:mRNA')){
					$validation_string .= "mRNA:$mRNA_ID|$validation_text;";
					$mRNA_validation_status += $validation_error_code{GFF_mRNA};
				}
				
				my @exons = $mRNA->get_SeqFeatures('exon');
				my $exon_error = 0;
				foreach my $exon (@exons){
					my %exon_hash;
					my %exon_attb = $exon->attributes;
					my $exon_ID = $exon_attb{load_id}->[0];
					$exon_hash{scaffold} = $exon->seq_id;
					$exon_hash{strand}   = $exon->strand;
					$exon_hash{start}    = $exon->start;
					$exon_hash{end}      = $exon->end;
				
					if(!_check_subfeature(\$validation_text,$exon_ID,\%mRNA_hash,\%exon_hash,'mRNA:exon')){
						$validation_string .=  "exon:$exon_ID|$validation_text;";
						$exon_error = 1;
					}
				
				}
				
				if($exon_error){ $mRNA_validation_status += $validation_error_code{GFF_exon}; } 
				
				my $CDS_sequence   = Initialize::get_CDS($mRNA,$validation_fh);
				unless($CDS_sequence){$errorLog->error("No CDS for $mRNA_ID");}
				my $prot_seq       = Initialize::get_translation($CDS_sequence,$mRNA_hash{strand});
				unless($prot_seq){$errorLog->error("No translation for $mRNA_ID");}
				my $sequence_check = _check_seq($prot_seq,\$validation_text);
				
				if($sequence_check ne 'passed'){										
					$validation_string .= "mRNA:$mRNA_ID";
					$validation_string .= "$validation_text;";
					$errorLog->error($validation_string);
					my($start_flag,$stop_flag,$internal_stop_count) = split /:/,$sequence_check;
					if(!$start_flag){
						$mRNA_validation_status += $validation_error_code{CDS_start};
					}
					
					if(!$stop_flag){
						$mRNA_validation_status += $validation_error_code{CDS_stop};
					}
					
					if(!$internal_stop_count == 0){
						$mRNA_validation_status += $validation_error_code{CDS_internal_stop};
					}
					
				}
				
				$mRNA->add_tag_value('validation_error_code' => $mRNA_validation_status);
				$mRNA->update();
		};
    #warn "$mRNA_validation_status : $validation_error_code{CDS_stop} : $NO_STOP";
		if($@){
			$errorLog->error($@);
			$gene_validation_status = -1;
		}elsif($mRNA_validation_status and $gene_validation_status > -1 and $gene_validation_status < $mRNA_validation_status){
			if($mRNA_validation_status == $validation_error_code{CDS_start} and  $NO_ATG == 2){
				#approved
				$gene_validation_status = 0;
			}elsif($mRNA_validation_status == $validation_error_code{CDS_stop} and $NO_STOP == 2){
				#approved
				$gene_validation_status = 0;
			}elsif(($mRNA_validation_status == ($validation_error_code{CDS_stop} + $validation_error_code{CDS_start})) and $NO_ATG == 2 and $NO_STOP == 2){
				#approved
				$gene_validation_status = 0;
			}elsif($mRNA_validation_status == $validation_error_code{CDS_start} and $NO_ATG == 1){
				$gene_validation_status = $mRNA_validation_status;				
			}elsif($mRNA_validation_status == $validation_error_code{CDS_stop} and $NO_STOP == 1){
				$gene_validation_status = $mRNA_validation_status;	
			}elsif(($mRNA_validation_status == ($validation_error_code{CDS_stop} + $validation_error_code{CDS_start})) and $NO_ATG == 1 and $NO_STOP == 1){
				$gene_validation_status = $mRNA_validation_status;
			}else{$gene_validation_status = -$mRNA_validation_status;}
		}		
	}
	
	if($gene_validation_status){
		my $gene_id = $gene_attb{load_id}->[0];
		my $owner =   $gene_attb{owner}->[0];
		
		if(!$owner){
			$owner = 'None';
		}
		
		print $validation_fh "gene:$gene_id\t$owner\t$validation_string\n";
	
	}
	$gene->add_tag_value('validation_error_code' => $gene_validation_status);
	return $gene_validation_status;
}

sub _check_subfeature{
	my($validation_text,$ID,$feature,$sub_feature,$key) = @_;
	my $passed = 1;	
	my ($feature_name,$sub_feature_name) = split /:/,$key;
	$$validation_text='';
	if($feature->{scaffold} ne $sub_feature->{scaffold}){
		$$validation_text .= "error:Scaffold not identical|proof:[$feature_name scaffold $feature->{scaffold}, $sub_feature_name scaffold $sub_feature->{scaffold}];";
		$passed = 0;
	}
	
	if($feature->{strand} ne $sub_feature->{strand}){
		$$validation_text .= "error:Strand not identical|proof:[$feature_name strand $feature->{strand}, $sub_feature_name strand $sub_feature->{strand}];";
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

sub _check_seq {
	my($prot_seq,$validation_text) = @_;
	my $start_aa = substr($prot_seq,0,1);
	my $end_aa   = substr($prot_seq,-1,1);
	my $start_flag=0;
	my $stop_flag =0;
	my $internal_stop_count =0;
	if($start_aa eq 'M'){$start_flag++;};
	if($end_aa eq '*'){$stop_flag++;};
	$internal_stop_count -= $stop_flag;
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
		return 'passed';
	}else{return "$start_flag:$stop_flag:$internal_stop_count";}
}
1;
 