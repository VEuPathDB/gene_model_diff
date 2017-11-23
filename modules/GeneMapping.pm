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

package GeneMapping;

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use TranscriptMapping;
use GeneModel;

sub resolve_maptype_cluster {
	my($dbh) = @_;
	
	my $clusters = GeneClusters::get_distinct_cluster_ids($dbh);
	
	foreach my $cluster (@{$clusters}){
		my($cluster_id) = @{$cluster};
		my($cap_gene_count,$cap_transcript_count,$vb_gene_count,$vb_transcript_count) = @{GeneClusters::get_cluster_summary_by_id($dbh,$cluster_id)->[0]}; 
		
		if($cap_gene_count == 1 and $vb_gene_count == 0){    #new Gene
			my($cap_gene) = @{GeneClusters::get_gene_cluster_by_id($dbh,$cluster_id)->[0]};
			insert_gene_mappings($dbh,"$cap_gene:",'new');
		}elsif($cap_gene_count == 1 and $vb_gene_count == 1){#gene change or identical gene
			my $cluster_ref = GeneClusters::get_gene_cluster_by_id($dbh,$cluster_id);
			my($cap_gene_id,$vb_gene_id,$final_maptype);
			
			foreach my $row (@{$cluster_ref}){
				my($gene_id,$source) = @{$row};
				
				if($source eq 'cap'){
					$cap_gene_id = $gene_id;
				}elsif($source eq 'vb'){
					$vb_gene_id = $gene_id;
				}else{croak("Source not known");}
			}
			
			my %transcript_lookup;
			my $transcript_maps = TranscriptMapping::get_all_transcript_mappings($dbh);
			foreach my $row (@{$transcript_maps}){
				my($cap_trans_id,$vb_trans_id,$map_type) = @{$row};
				$transcript_lookup{$cap_trans_id} = $map_type;
				if(defined($vb_trans_id)){
					$transcript_lookup{$vb_trans_id}  = $map_type;
				}
			}
			
			
			my @transcript = (@{GeneModel::get_gene_model_by_id($dbh,'gene',$cap_gene_id,'cap','transcript',1)},@{GeneModel::get_gene_model_by_id($dbh,'gene',$vb_gene_id,'vb','transcript',1)});
			
			
			foreach my $transcript (@transcript){			
				my $maptype = $transcript_lookup{$transcript};
				if(!defined($maptype) or $maptype ne 'identical'){
					$final_maptype = 'change';
				}else{$final_maptype = 'identical';}
			}
			
			insert_gene_mappings($dbh,"$cap_gene_id:$vb_gene_id",$final_maptype);
		
		}elsif($cap_gene_count > 1 and $vb_gene_count == 1){#split
			my $cluster_ref = GeneClusters::get_gene_cluster_by_id($dbh,$cluster_id);
			my $vb_gene_id;
			my @cap_gene_id;
			foreach my $row (@{$cluster_ref}){
				my ($gene_id,$source) = @{$row};
				if($source eq 'vb'){
					$vb_gene_id = $gene_id;
				}else{
					push @cap_gene_id,$gene_id; 	
				}
			}
			
			foreach my $cap_gene_id (@cap_gene_id){
				insert_gene_mappings($dbh,"$cap_gene_id:$vb_gene_id",'split');
			}
			
		}elsif($cap_gene_count == 1 and $vb_gene_count > 1){#merge
			my $cluster_ref = GeneClusters::get_gene_cluster_by_id($dbh,$cluster_id);
			my $cap_gene_id;
			my @vb_gene_id;
			foreach my $row (@{$cluster_ref}){
				my ($gene_id,$source) = @{$row};
				if($source eq 'cap'){
					$cap_gene_id = $gene_id;
				}else{
					push @vb_gene_id,$gene_id; 	
				}
			}
			
			foreach my $vb_gene_id (@vb_gene_id){
				insert_gene_mappings($dbh,"$cap_gene_id:$vb_gene_id",'merge');
			}
		
		}elsif($cap_gene_count > 1 and $vb_gene_count > 1){ #complex loci i.e merge and split
			my $cluster_ref = GeneClusters::get_gene_cluster_by_id($dbh,$cluster_id);
			my %gene_lookup;
			
			my $transcript_maps = TranscriptMapping::get_all_transcript_mappings($dbh);			
			foreach my $row (@{$transcript_maps}){
				my($cap_trans_id,$vb_trans_id,$map_type) = @{$row};
				next unless($cap_trans_id and $vb_trans_id);
				my $cap_gene_id = GeneModel::get_distinct_id_by_id($dbh,'transcript',$cap_trans_id,'gene','cap',0);
				my $vb_gene_id  = GeneModel::get_distinct_id_by_id($dbh,'transcript',$vb_trans_id,'gene','vb',0);
				$gene_lookup{cap}{$cap_gene_id}{count}++;
				$gene_lookup{vb}{$vb_gene_id}{count}++;
				push @{$gene_lookup{vb}{$vb_gene_id}{genes}},$cap_gene_id;
				push @{$gene_lookup{cap}{$cap_gene_id}{genes}},$vb_gene_id;
			}
			
			my %uniq_gene_mappings;
			foreach my $row (@{$cluster_ref}){
				my($gene_id,$source) = @{$row};				
				if($source eq 'vb'){										
					if($gene_lookup{vb}{$gene_id}{count} > 1){
						foreach my $cap_gene_id (@{$gene_lookup{vb}{$gene_id}{genes}}){
							if(! exists $uniq_gene_mappings{split}{"$cap_gene_id:$gene_id"}){
								insert_gene_mappings($dbh,"$cap_gene_id:$gene_id",'split');
								$uniq_gene_mappings{split}{"$cap_gene_id:$gene_id"} = 1;
							}
						}
					}
				}	
			
				
				if($source eq 'cap'){
					if($gene_lookup{cap}{$gene_id}{count} > 1){
						foreach my $vb_gene_id (@{$gene_lookup{cap}{$gene_id}{genes}}){
							if(! exists $uniq_gene_mappings{merge}{"$gene_id:$vb_gene_id"}){
								insert_gene_mappings($dbh,"$gene_id:$vb_gene_id",'merge');
								$uniq_gene_mappings{merge}{"$gene_id:$vb_gene_id"} = 1;
							}
						}
					}						
				}				
			 }
		
		}else{confess("Cluster $cluster_id was not processed");}#error
	}
		
}

sub insert_gene_mappings{
	my($dbh,$gene_ids,$map_type) = @_;
	my $insert_gene_mappings_sth = $dbh->prepare(get_sql('insert_gene_mappings'));
	my($cap_gene_id,$vb_gene_id) = split /:/,$gene_ids;
	
	if($map_type eq 'new'){
		$vb_gene_id = undef;			
	}	
	$insert_gene_mappings_sth->execute($cap_gene_id,$vb_gene_id,$map_type);
}



sub _add_to_hash {
	my($hash_ref,$array_ref,$colunm) = @_;
	foreach my $row (@{$array_ref}){
		if(!exists $hash_ref->{$row->[$colunm]}){
			$hash_ref->{$row->[$colunm]} = 1;
		}	
	}
	
}

sub update_maptype {
	my ($dbh,$id_type,$id,$maptype) = @_;
	
	my $sql = "update gene_mappings set map_type = \'$maptype\' where $id_type = \'$id\';"; 
	$dbh->do($sql);

}

sub get_gene_mappings_by_maptype {
	my ($dbh,$maptype) = @_;
	
	my $sql = "select cap_gene_id,vb_gene_id from gene_mappings where map_type = \'$maptype\';";
	my $array_ref = _submit_sql($dbh,$sql);	
	return $array_ref;
	
}

sub get_gene_mappings_by_id {
	my ($dbh,$type,$id,$return_id) = @_;
	
	my $sql = "select $return_id from gene_mappings where $type = \'$id\';";
	my $array_ref = _submit_sql($dbh,$sql);	
	return $array_ref;	
}

sub get_gene_mappings_by_id_and_maptype {
	my ($dbh,$type,$id,$return_id,$maptype) = @_;
	
	my $sql = "select $return_id from gene_mappings where $type = \'$id\' and map_type = \'$maptype\';";
	my $array_ref = _submit_sql($dbh,$sql);	
	return $array_ref;	
}

sub get_gene_mappings_by_group {
	my ($dbh,$gene_id_A,$gene_id_B,$map_type) = @_;
    
	my $sql = "select $gene_id_A from gene_mappings where map_type = \'$map_type\' group by $gene_id_A having count(distinct $gene_id_B) > 1;";
	my $array_ref = _submit_sql($dbh,$sql);
	return $array_ref;
}

sub _submit_sql {
	my ($dbh,$sql) = @_;
	my $array_ref = $dbh->selectall_arrayref($sql);
	return $array_ref;	
}

sub get_sql{
	my ($sql_name) = @_;

	if($sql_name eq 'insert_gene_mappings'){
			my $insert_gene_mappings = "insert gene_mappings(cap_gene_id,vb_gene_id,map_type) select ?,?,?;";
			return $insert_gene_mappings;
	}
}
1;