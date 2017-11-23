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

package TranscriptLinks;

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use GeneModel;
use lib('.');
use ExonMapping;

sub work_out_transcript_links {
	my($dbh) = @_;
	
	my $gene_clusters = GeneClusters::get_distinct_cluster_ids($dbh);
	my $insert_transcript_links_sth = $dbh->prepare(get_sql('insert_transcript_links'));
	foreach my $cluster_row (@{$gene_clusters}){
		my($cluster_id) = @{$cluster_row};
		my $cluster_summary = GeneClusters::get_cluster_summary_by_id($dbh,$cluster_id);
		my($cap_gene_count,$cap_transcript_count,$vb_gene_count,$vb_transcript_count) = @{$cluster_summary->[0]};
		
		if($cap_gene_count == 1 and $vb_gene_count == 1 and ($cap_transcript_count > 1 or $vb_transcript_count > 1)){
				  link_matched_transcripts($dbh,$cluster_id);
		}elsif($cap_gene_count > 0 and $vb_gene_count > 0){ 
				link_overlapping_transcripts($dbh,$cluster_id);		
		}elsif($cap_gene_count == 1 and $vb_gene_count == 0){ #New Gene
				my $cluster_ref = GeneClusters::get_gene_cluster_by_id($dbh,$cluster_id);
				my($gene_id,$source) = @{$cluster_ref->[0]};
				my $transcript_ref = GeneModel::get_distinct_id_by_id($dbh,'gene',$gene_id,'transcript',$source,1);
				foreach my $row (@{$transcript_ref}){
					my ($transcript_id) = @{$row};
					$insert_transcript_links_sth->execute($transcript_id,undef);
				}
		}else{  #Error 
				croak("Cluster must atleast have one cap gene");
		}
	
	}
	

}

sub link_matched_transcripts {
	my($dbh,$cluster_id) = @_;	
	my $insert_transcript_links_sth = $dbh->prepare(get_sql('insert_transcript_links'));
	my $cluster_ref = GeneClusters::get_gene_cluster_by_id($dbh,$cluster_id);
	my %linked_transcripts;
	my %allready_linked;
	my %unmapped_transcripts;
	foreach my $row (@{$cluster_ref}){		
		my($gene_id,$source) = @{$row};
		my $transcript_ref = GeneModel::get_distinct_id_by_id($dbh,'gene',$gene_id,'transcript',$source,1);
		foreach my $row (@{$transcript_ref}){
			my ($transcript_id) = @{$row};
			my $transcript_exon_mappings_array_ref = ExonMapping::get_exon_mappings_by_id($dbh,$source,$transcript_id);
			foreach my $row (@{$transcript_exon_mappings_array_ref}){
				my($cap_exon_id,$vb_exon_id,$map_type) = @{$row};
				next unless defined($vb_exon_id);
				my $cap_transcript_id = GeneModel::get_gene_model_by_id($dbh,'exon',$cap_exon_id,'cap','transcript');
				my $vb_transcript_id  = GeneModel::get_gene_model_by_id($dbh,'exon',$vb_exon_id,'vb','transcript');
				my $link_key = "$cap_transcript_id:$vb_transcript_id";
				
				if($map_type eq 'identical'){
					$linked_transcripts{$link_key}++;
				}								
			}
		}
	}
	
	#transcripts is linked to eachother based on most identical exons. 
	foreach my $link_key (sort {$linked_transcripts{$b} <=> $linked_transcripts{$a}} keys %linked_transcripts ){
		
		my($cap_transcript_id,$vb_transcript_id) = split /:/,$link_key;
		if(!exists $allready_linked{$cap_transcript_id} and !exists $allready_linked{$vb_transcript_id}){
			$allready_linked{$cap_transcript_id} = 1;
			$allready_linked{$vb_transcript_id}  = 1;
			$insert_transcript_links_sth->execute($cap_transcript_id,$vb_transcript_id);
		}elsif(!exists $unmapped_transcripts{$link_key}){
			$unmapped_transcripts{$link_key} = $linked_transcripts{$link_key};
		}
	}
	
	foreach my $link_key (keys %unmapped_transcripts){
		my($cap_transcript_id,$vb_transcript_id) = split /:/,$link_key;
		if(!exists $allready_linked{$cap_transcript_id} or !exists $allready_linked{$vb_transcript_id}){
			$allready_linked{$cap_transcript_id} = 1;
			$allready_linked{$vb_transcript_id}  = 1;
			$insert_transcript_links_sth->execute($cap_transcript_id,$vb_transcript_id);
		}
	
	}	
}

sub link_overlapping_transcripts {
	my($dbh,$cluster_id) = @_;
	my $insert_transcript_links_sth = $dbh->prepare(get_sql('insert_transcript_links'));
	my $cluster_ref = GeneClusters::get_gene_cluster_by_id($dbh,$cluster_id);
	my %linked_transcripts;
	foreach my $row (@{$cluster_ref}){
		
		my($gene_id,$source) = @{$row};
		my $transcript_ref = GeneModel::get_distinct_id_by_id($dbh,'gene',$gene_id,'transcript',$source,1);
		foreach my $row (@{$transcript_ref}){
			my ($transcript_id) = @{$row};
			my $transcript_exon_mappings_array_ref = ExonMapping::get_exon_mappings_by_id($dbh,$source,$transcript_id);
			foreach $row (@{$transcript_exon_mappings_array_ref}){
				my($cap_exon_id,$vb_exon_id,$map_type) = @{$row};
				next unless(defined($cap_exon_id) and defined($vb_exon_id));
				if($cap_exon_id eq 'none' or $vb_exon_id eq 'none'){next;} 
				my $cap_transcript_id= GeneModel::get_gene_model_by_id($dbh,'exon',$cap_exon_id,'cap','transcript');
				my $vb_transcript_id = GeneModel::get_gene_model_by_id($dbh,'exon',$vb_exon_id,'vb','transcript');
				if(! exists $linked_transcripts{$cap_transcript_id}{$vb_transcript_id} ){
								$linked_transcripts{$cap_transcript_id}{$vb_transcript_id} = $map_type;
								$insert_transcript_links_sth->execute($cap_transcript_id,$vb_transcript_id);
				}
			}
		}
	}
}

sub get_all_linked_transcripts{
	my($dbh) =@_;
	my $sql = "select cap_transcript_id,vb_transcript_id from transcript_links;";
	my $array_ref = $dbh->selectall_arrayref($sql);
	return $array_ref;		
}

sub get_sql{
	my ($sql_name) = @_;
		
	if($sql_name eq 'insert_transcript_links'){
		my $insert_transcript_links = "insert transcript_links(cap_transcript_id,vb_transcript_id) select ?,?";
		return $insert_transcript_links;
	}elsif($sql_name eq 'get_maped_unlinked_vb_transcripts'){
		my $maped_unlinked_vb_transcripts = "select distinct transcript_id from gene_model,exon_mappings 
											 where exon_id =  vb_exon_id
											 and source = 'vb'
											 and not exists (select * from transcript_links where vb_transcript_id = transcript_id);";
	    return $maped_unlinked_vb_transcripts;
	}	
}
1;