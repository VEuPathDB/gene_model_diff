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
package AnnotationEvents;

use strict;
use warnings;

sub get_identical_gene {
	my($dbh) = @_;
	
	my $gene_count=0,
	my $identical_gene_insert_sth = $dbh->prepare(get_sql('gene_events'));
	
	
	my $identical_gene_aref = GeneMapping::get_gene_mappings_by_maptype($dbh,'identical'); 
	
	foreach my $row (@{$identical_gene_aref}){
			my $vb_gene_id  = $row->[1];
		    my $cap_gene_id = $row->[0];
			$identical_gene_insert_sth->execute($vb_gene_id,$cap_gene_id,'identical');
			$gene_count++;
	}		
	
	return $gene_count;
}

sub get_new_gene{
	my($dbh) = @_;
	
	my $gene_count=0;
	my $new_gene_insert_sth = $dbh->prepare(get_sql('gene_events'));
	
	
	my $new_gene_aref = GeneMapping::get_gene_mappings_by_maptype($dbh,'new');
	
	foreach my $row (@{$new_gene_aref}){
			my $cap_gene_id = $row->[0];
			$new_gene_insert_sth->execute(undef,$cap_gene_id,'new_gene');
			$gene_count++;
	}
	return $gene_count;
}

sub get_changed_genes {
	my($dbh) = @_;
	
	my $gene_count=0;
	my $changed_gene_insert_sth = $dbh->prepare(get_sql('gene_events'));
	
	
	my $changed_gene_aref = GeneMapping::get_gene_mappings_by_maptype($dbh,'change');
	
	foreach my $row (@{$changed_gene_aref}){
			my $vb_gene_id  = $row->[1];
		    my $cap_gene_id = $row->[0];
			$changed_gene_insert_sth->execute($vb_gene_id,$cap_gene_id,'change_gene');
			$gene_count++;
	}
	
	return $gene_count;
}

sub get_splits{
	my($dbh) = @_;
	
	my($total_splits,$vb_genes,$cap_genes);
	$total_splits=$vb_genes=$cap_genes=0;
	my $events_insert_sth = $dbh->prepare(get_sql('gene_events'));
	
	
	my $splits_array_ref = GeneMapping::get_gene_mappings_by_group($dbh,'vb_gene_id','cap_gene_id','split');
	
	foreach my $row (@{$splits_array_ref}){
			my $vb_gene_id = $row->[0];
			$total_splits++;
			$vb_genes--;
			my $cap_id_string;
			my $cap_gene_array = GeneMapping::get_gene_mappings_by_id_and_maptype($dbh,'vb_gene_id',$vb_gene_id,'cap_gene_id','split');
			foreach my $row (@{$cap_gene_array}){
					my $cap_gene_id = $row->[0];
					$cap_id_string .= $cap_gene_id . ':';
					$cap_genes++;
			}
			$cap_id_string=~s/\:$//;
			$events_insert_sth->execute($vb_gene_id,$cap_id_string,'split_gene');
	}
	
	return ($total_splits,$vb_genes,$cap_genes);
}

sub get_merge{
	my($dbh) = @_;
	
	my ($total_merge,$vb_genes,$cap_genes);
	$total_merge=$vb_genes=$cap_genes=0;
	my $events_insert_sth = $dbh->prepare(get_sql('gene_events'));
	
	my $cap_merge_aref = GeneMapping::get_gene_mappings_by_group($dbh,'cap_gene_id','vb_gene_id','merge');
	
	foreach my $row (@{$cap_merge_aref}){
		my $cap_gene_id = $row->[0];
		
		$total_merge++;
		$cap_genes++;
		
		my $vb_id_string;
		my $vb_merge_aref = GeneMapping::get_gene_mappings_by_id_and_maptype($dbh,'cap_gene_id',$cap_gene_id,'vb_gene_id','merge');
		foreach my $row (@{$vb_merge_aref}){
			my $vb_gene_id = $row->[0];
			$vb_id_string .= $vb_gene_id . ':';
			$vb_genes--;
		}
		$vb_id_string=~s/\:$//;
		$events_insert_sth->execute($vb_id_string,$cap_gene_id,'merge_gene');
	}
	
	return ($total_merge,$vb_genes,$cap_genes);	
}

sub get_sql{
	my ($sql_name) = @_;
	if($sql_name eq 'gene_events'){
		    my $gene_events = "insert gene_events(vb_gene_id,cap_gene_id,events) select ?,?,?";
		    return $gene_events;
	}
}
1;