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

package TranscriptMapping;
use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use TranscriptLinks;
use ExonMapping;
use CDS;

sub resolve_transcript_mappings{
	my($dbh) = @_;		
	my $insert_transcript_mappings_sth = $dbh->prepare(get_sql('insert_transcript_mappings'));
		
	my $linked_transcript_array_ref = TranscriptLinks::get_all_linked_transcripts($dbh);
	foreach my $linked_pair (@{$linked_transcript_array_ref}){
				my ($cap_transcript_id,$vb_transcript_id) = @{$linked_pair};
				my %linked_transcripts;
				#warn "$cap_transcript_id , $vb_transcript_id";
				
				$linked_transcripts{cap}{id} = $cap_transcript_id;
				
				my $cap_exons_for_transcript_array_ref = GeneModel::get_gene_model_by_id($dbh,'transcript',$cap_transcript_id,'cap');
				
				foreach my $row (@{$cap_exons_for_transcript_array_ref}){					
					$linked_transcripts{cap}{exons}{$row->[0]}=1;	
				}
				my $vb_exons_for_transcript_array_ref;				
				if(defined($vb_transcript_id)){
					$vb_exons_for_transcript_array_ref = GeneModel::get_gene_model_by_id($dbh,'transcript',$vb_transcript_id,'vb');
				
					foreach my $row (@{$vb_exons_for_transcript_array_ref}){					
						$linked_transcripts{vb}{exons}{$row->[0]}=1;
					}
				}
				
				my $map_type = compare_exons($dbh,\%linked_transcripts);
				#warn "MAP_TYPE: $map_type";
				#If all exons are identical we check cds
				
				if($map_type eq 'identical'){
					my ($cap_cds_start,$cap_cds_end) = CDS::get_cds_pos_by_parent_id($dbh,$cap_transcript_id);
					my ($vb_cds_start,$vb_cds_end) = CDS::get_cds_pos_by_parent_id($dbh,$vb_transcript_id);
					if($cap_cds_start ne $vb_cds_start or $cap_cds_end ne $vb_cds_end){
						$map_type = 'change';
					}
				}
				$insert_transcript_mappings_sth->execute($cap_transcript_id,$vb_transcript_id,$map_type);
	}	
}

sub compare_exons {
	my($dbh,$linked_transcripts_ref) = @_;
	my $errorLog = Log::Log4perl->get_logger("TranscriptMapping::compare_exons");
	my $resulting_cap_map_type;
	my $resulting_vb_map_type;
	my $cap_mapped_to_exon_in_other_transcript;
	my $vb_mapped_to_exon_in_other_transcript;
	foreach my $exon_id (keys %{$linked_transcripts_ref->{cap}{exons}}){
		#warn "CAP EXON ID $exon_id";
		my $EXON_map_type=undef;
		
		my $cap_exon_mappings_array_ref = ExonMapping::get_exon_mappings($dbh,'cap',$exon_id);
		foreach my $row (@{$cap_exon_mappings_array_ref}){
			my($map_type,$vb_exon_id) = @{$row};
			#warn "comparing cap:$exon_id vb:$vb_exon_id Type:$map_type ";
			#$EXON_map_type = $map_type;
			if(!$map_type or !$vb_exon_id or ($vb_exon_id eq 'none')){
				$resulting_cap_map_type = 'new';
				$EXON_map_type = 'new';
			}elsif(exists $linked_transcripts_ref->{vb}{exons}{$vb_exon_id}){
				$EXON_map_type = $map_type;
				if(!$resulting_cap_map_type){ 
					$resulting_cap_map_type = $map_type;
					$EXON_map_type  = $map_type;
				}elsif($resulting_cap_map_type ne $map_type){
					$resulting_cap_map_type = 'change';
					$EXON_map_type = 'change';
				}
			}else{  #warn "This Exon does not exists in transcript";
					$cap_mapped_to_exon_in_other_transcript = 1; }			
					#warn "EXON_map_type 1 $EXON_map_type";
		}
		
		#warn "EXON_map_type 2 $EXON_map_type";
		if((!$EXON_map_type or !$resulting_cap_map_type) and $cap_mapped_to_exon_in_other_transcript){
			$resulting_cap_map_type = 'change';		
		}elsif(!$resulting_cap_map_type){$errorLog->error("CAP exon: $exon_id was linked, but the maptype >$EXON_map_type< could not be resolved.");}
		
		#warn "comparing cap to vb: $resulting_cap_map_type";
	}
	
	
	
	foreach my $exon_id (keys %{$linked_transcripts_ref->{vb}{exons}}){
		#warn "VB EXON ID $exon_id";
		my $EXON_map_type=undef;
		  
		my $vb_exon_mappings_array_ref = ExonMapping::get_exon_mappings($dbh,'vb',$exon_id);
		foreach my $row (@{$vb_exon_mappings_array_ref}){
			my($map_type,$cap_exon_id) = @{$row};
			#warn "comparing vb:$exon_id cap:$cap_exon_id Type:$map_type ";
			#$EXON_map_type = $map_type;
			if(!$map_type or !$cap_exon_id or ($cap_exon_id eq 'none')){
				$resulting_cap_map_type = 'change';
				$EXON_map_type = 'change';
			}elsif(exists $linked_transcripts_ref->{cap}{exons}{$cap_exon_id}){
				$EXON_map_type = $map_type;
				if(!$resulting_vb_map_type){ 
					$resulting_vb_map_type = $map_type;
					$EXON_map_type = $map_type;
				}elsif($resulting_vb_map_type ne $map_type){
					$resulting_vb_map_type = 'change';
					$EXON_map_type = 'change';
				}
			}else{ 	#warn "This Exon does not exists in transcript";
					$vb_mapped_to_exon_in_other_transcript = 1;}
		}
		
		if((!$EXON_map_type or !$resulting_vb_map_type) and $vb_mapped_to_exon_in_other_transcript){
			$resulting_vb_map_type = 'change';		
		}elsif(!$resulting_vb_map_type){$errorLog->error("VB exon: $exon_id was linked, but the maptype >$EXON_map_type< could not be resolved.");}
		
		#warn "comparing vb to cap: $resulting_vb_map_type";
	}
	

	
	#warn "RESULT $resulting_cap_map_type $resulting_vb_map_type";
	
	if(!defined($resulting_vb_map_type) or ($resulting_cap_map_type eq $resulting_vb_map_type)){
		return $resulting_cap_map_type;	
	}elsif(($resulting_cap_map_type eq $resulting_vb_map_type) or !$resulting_cap_map_type){
		return $resulting_vb_map_type;
	}else{
		return 'change';	
	}
	
}

sub get_all_transcript_mappings {
	my ($dbh) = @_;
	my $sql = "select cap_trans_id,vb_trans_id,map_type from transcript_mappings;";
	my $array_ref = $dbh->selectall_arrayref($sql);
	return $array_ref;
}

sub insert_transcript_mappings{
		my($dbh)=@_;
		
	    
		my $insert_transcript_mappings_sth = $dbh->prepare(get_sql('insert_transcript_mappings'));
		
		#no_match_cap
		
		
		
		my $exon_array_ref = ExonMapping::get_exon_mappings($dbh);
		my %transcript_mapping;
		my %seen_transcript;
		foreach my $row (@{$exon_array_ref}){
			my($cap_exon_id,$vb_exon_id,$map_type) = ($row->[0],$row->[1],$row->[2]);
			
			
			
			
			my $cap_transcript_id = GeneModel::get_gene_model_by_id($dbh,'exon',$cap_exon_id,'cap','transcript');
			my $vb_transcript_id  = GeneModel::get_gene_model_by_id($dbh,'exon',$cap_exon_id,'cap','transcript');
			#warn "CAP $cap_transcript_id - VB $vb_transcript_id\n";
			if(!$vb_transcript_id){
				$vb_transcript_id ='none';
			}
			if(!$cap_transcript_id){
				foreach my $cap_transcript_key (keys %transcript_mapping){
					foreach my $vb_transcript_value (@{$transcript_mapping{$cap_transcript_key}{vb_transcript_id}}){
							if($vb_transcript_value eq $vb_transcript_id){
								$transcript_mapping{$cap_transcript_key}{macth_type} = 'change';
							}
					}
				}
				next;				
			}
			#my $key = $cap_transcript_id . ':' . $vb_transcript_id;
			
			
			my $macth_type;
			if($row->[2] eq 'no_match_cap'){
				$macth_type= 'new';
			}elsif($row->[2] eq 'no_match_vb'){
				$macth_type= 'lost';
			}elsif($row->[2] eq 'identical'){
				$macth_type= 'identical';
			}else{
				$macth_type= 'change';	
			}
			
			
			
			
			if(!exists $transcript_mapping{$cap_transcript_id}){
				$transcript_mapping{$cap_transcript_id}{macth_type} = $macth_type;
				if($vb_transcript_id ne 'none' and !exists $seen_transcript{$cap_transcript_id}{$vb_transcript_id}){
						push @{$transcript_mapping{$cap_transcript_id}{vb_transcript_id}},$vb_transcript_id;
						$seen_transcript{$vb_transcript_id} = 1;				
				}
			}else{
				
				if($transcript_mapping{$cap_transcript_id}{macth_type} ne $macth_type){
					$transcript_mapping{$cap_transcript_id}{macth_type} = 'change';
				}
				
				foreach my $transcript_id (@{$transcript_mapping{$cap_transcript_id}{vb_transcript_id}}){
						if($transcript_id ne $vb_transcript_id){
							$transcript_mapping{$cap_transcript_id}{macth_type} = 'change';
							if($vb_transcript_id ne 'none' and !exists $seen_transcript{$cap_transcript_id}{$vb_transcript_id}){
								push @{$transcript_mapping{$cap_transcript_id}{vb_transcript_id}},$vb_transcript_id;
								$seen_transcript{$vb_transcript_id} = 1;
							}
						}
				}	
			}
			
		}
		
		foreach my $cap_transcript_id (keys %transcript_mapping){
			if(!defined(@{$transcript_mapping{$cap_transcript_id}{vb_transcript_id}}) or scalar @{$transcript_mapping{$cap_transcript_id}{vb_transcript_id}} == 0){
				$insert_transcript_mappings_sth->execute($cap_transcript_id,'none',$transcript_mapping{$cap_transcript_id}{macth_type});		
			}
			warn scalar @{$transcript_mapping{$cap_transcript_id}{vb_transcript_id}};
			if(scalar @{$transcript_mapping{$cap_transcript_id}{vb_transcript_id}} > 1 ){
				$transcript_mapping{$cap_transcript_id}{macth_type} = 'change';		
			}
			foreach my $vb_transcript_id (@{$transcript_mapping{$cap_transcript_id}{vb_transcript_id}}){
				$insert_transcript_mappings_sth->execute($cap_transcript_id,$vb_transcript_id,$transcript_mapping{$cap_transcript_id}{macth_type});
			}
		}
		
}

sub get_sql{
	my ($sql_name) = @_;
	
	if($sql_name eq 'insert_transcript_mappings'){
			my $insert_transcript_mappings = "insert transcript_mappings(cap_trans_id,vb_trans_id,map_type) select ?,?,?;";
			return $insert_transcript_mappings;
	}
}
1;