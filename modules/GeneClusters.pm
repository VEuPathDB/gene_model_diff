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

package GeneClusters;

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use GeneModel;
use lib('.');
use ExonMapping;

sub work_out_gene_clusters{
	my($dbh) = @_;
	my %added_cap_genes;
	my $cap_gene_model_ref = GeneModel::get_distinct_id_by_source($dbh,'gene','cap');
	foreach my $cap_gene_id_ref (@{$cap_gene_model_ref}){
		if(!exists $added_cap_genes{$cap_gene_id_ref->[0]}){
			my %gene_cluster;			
			recursive_cluster_genes($dbh,$cap_gene_id_ref,\%gene_cluster);
			insert_gene_cluster($dbh,\%gene_cluster);
			foreach my $cap_gene_id (keys $gene_cluster{cap}){				
				$added_cap_genes{$cap_gene_id} = 1; 
			}		
		}		
	}	
}

sub calculate_cluster_summary {
	my($dbh) = @_;
	my $clusters_ref = get_gene_clusters($dbh);
	my %clusters;
	foreach my $cluster_row (@{$clusters_ref}){
		my($cluster_id,$gene_id,$source) = @{$cluster_row};
		
		my $transcripts = GeneModel::get_distinct_id_by_id($dbh,'gene',$gene_id,'transcript',$source,1);
		
		$clusters{$cluster_id}{$source}{$gene_id} = scalar(@{$transcripts});				
	}
	
	foreach my $cluster_id (keys %clusters){
		my($cap_gene_count,$cap_trans_count,$vb_gene_count,$vb_trans_count);		
		foreach my $source (keys %{$clusters{$cluster_id}}){
			my $gene_count;
			my $transcript_count;
			foreach my $gene_id (keys %{$clusters{$cluster_id}{$source}}){
					$gene_count++;
					$transcript_count += $clusters{$cluster_id}{$source}{$gene_id};
			}
			if($source eq 'cap'){
				$cap_gene_count  = $gene_count;
				$cap_trans_count = $transcript_count;
			}
			elsif($source eq 'vb'){
				$vb_gene_count  = $gene_count;
				$vb_trans_count = $transcript_count;
			}
		}
	
		insert_cluster_summary($dbh,$cluster_id,$cap_gene_count,$cap_trans_count,$vb_gene_count,$vb_trans_count);
	}
}

sub recursive_cluster_genes{
    my($dbh,$cap_gene_array,$gene_cluster) = @_;
    if(scalar @{$cap_gene_array} == 0){
    	return 1;
    }
    
    my $cap_gene_id  = shift @{$cap_gene_array};
    if(!exists $gene_cluster->{cap}{$cap_gene_id}){
    	$gene_cluster->{cap}{$cap_gene_id} =1;
		my $vb_gene_ids  = get_gene_via_mapped_exons($dbh,$cap_gene_id,'cap','vb');
		foreach my $vb_gene_id (@{$vb_gene_ids}){
			if(!exists $gene_cluster->{vb}{$vb_gene_id}){
				$gene_cluster->{vb}{$vb_gene_id} = 1;
				my $cap_gene_ids = get_gene_via_mapped_exons($dbh,$vb_gene_id,'vb','cap');				
				push(@{$cap_gene_array},@{$cap_gene_ids});
			}
		}
    }
    recursive_cluster_genes($dbh,$cap_gene_array,$gene_cluster);
}

sub get_gene_via_mapped_exons{
	my($dbh,$query_gene_id,$query_source,$return_source) = @_;
	my %gene_id_duplicate;
	my @return_gene_ids;
	my @query_exon_ids = @{GeneModel::get_gene_model_by_id($dbh,'gene',$query_gene_id,$query_source,'exon',1)};
	foreach my $exon_id (@query_exon_ids){
			my $return_exon_ids = ExonMapping::get_exon_mappings($dbh,$query_source,$exon_id);
			foreach my $row (@{$return_exon_ids}){
				my $return_exon_id = @{$row}[1];
				next unless defined($return_exon_id);
				my $return_gene_id = GeneModel::get_distinct_id_by_id($dbh,'exon',$return_exon_id,'gene',$return_source);
				if(defined($return_gene_id) and !exists $gene_id_duplicate{$return_gene_id} ){
    	    		$gene_id_duplicate{$return_gene_id} =1;
					push @return_gene_ids,$return_gene_id;
    	    	}
			}	
	}
	return \@return_gene_ids;
}

sub get_gene_clusters{
	my($dbh) = @_;
	my $sql = "select gene_cluster_id,gene_id,source from gene_clusters;";
	my $array_ref = _submit_sql($dbh,$sql);
	return $array_ref;	
}

sub get_distinct_cluster_ids{
	my($dbh) = @_;
	my $sql = "select distinct(gene_cluster_id) from gene_clusters;";
	my $array_ref = _submit_sql($dbh,$sql);
	return $array_ref;
}

sub get_gene_cluster_by_id {
	my($dbh,$gene_cluster_id) = @_;
	my $sql = "select gene_id,source from gene_clusters where gene_cluster_id = $gene_cluster_id;";
	my $array_ref = _submit_sql($dbh,$sql);
	return $array_ref;
}

sub get_cluster_summary{
	my($dbh) = @_;
	my $sql = "select  gene_cluster_id,
			   cap_gene_count,
			   cap_transcript_count,
			   vb_gene_count,
			   vb_transcript_count
			   from cluster_summary;";
	my $array_ref = _submit_sql($dbh,$sql);		   
	return $array_ref;
}

sub get_cluster_summary_by_id{
	my($dbh,$cluster_id) = @_;
	my $sql = "select cap_gene_count,
			   cap_transcript_count,
			   vb_gene_count,
			   vb_transcript_count
			   from cluster_summary
			   where gene_cluster_id = $cluster_id;";
	my $array_ref = _submit_sql($dbh,$sql);		   
	return $array_ref;
}

sub insert_gene_cluster{
	my($dbh,$gene_cluster) = @_;	
	my $insert_sql = "insert gene_clusters(
			   gene_cluster_id, 
			   gene_id,
			   source)
			   select ?,?,?;";
	my $insert_sth = $dbh->prepare($insert_sql);
	#get cluster ID!
	my $get_last_cluster_id_sql = "select max(gene_cluster_id) from gene_clusters;";
	my $last_cluster_ref = _submit_sql($dbh,$get_last_cluster_id_sql);	
	my $last_cluster_id = defined($last_cluster_ref->[0]->[0])? $last_cluster_ref->[0]->[0] : 0;
	my $gene_cluster_id = $last_cluster_id + 1;
	foreach my $gene_id (keys %{$gene_cluster->{cap}}){		
		$insert_sth->execute($gene_cluster_id,$gene_id,'cap');
	}
	
	foreach my $gene_id (keys %{$gene_cluster->{vb}}){
		$insert_sth->execute($gene_cluster_id,$gene_id,'vb');
	}
	
}

sub insert_cluster_summary {
	my($dbh,$gene_cluster_id,$cap_gene_count,$cap_trans_count,$vb_gene_count,$vb_trans_count) = @_;
	
	if(!defined($vb_gene_count)){
		$vb_gene_count = 0;
		$vb_trans_count = 0;	
	}
	my $sql = "insert cluster_summary(
			   gene_cluster_id,
			   cap_gene_count,
			   cap_transcript_count,
			   vb_gene_count,
			   vb_transcript_count)
			   select ?,?,?,?,?";
			   
	my $insert_sth = $dbh->prepare($sql);
	$insert_sth->execute($gene_cluster_id,$cap_gene_count,$cap_trans_count,$vb_gene_count,$vb_trans_count);
}

sub _submit_sql {
	my($dbh,$sql) = @_;
	my $array_ref = $dbh->selectall_arrayref($sql);
	return $array_ref;	
}

1;