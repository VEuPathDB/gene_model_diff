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

AnnotationEvents

=head1 SYNOPSIS
	
	use AnnotationEvents;
	
	$identical_gene_count = AnnotationEvents::get_identical_gene($dbh);

=head1 DESCRIPTION

This module is the interface to the gene_events table which contains the final event type for each gene.
The module contains the SQL to insert the final event type based on the event type in the gene_mapping table. 

=head1 Author

	Mikkel B Christensen

=head1 METHODS

=cut

package AnnotationEvents;

use strict;
use warnings;

=head2 

 Title:	get_identical_gene    	
 Usage: AnnotationEvents::get_identical_gene($dbh)  	
 Function: Gets identical genes from the gene_mapping table and inserts them in to the gene_events table  	
 Returns: Count of indentical genes  
 Args: Database handle object     
=cut

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

=head2 

 Title: get_new_gene   	
 Usage: AnnotationEvents::get_new_gene($dbh)   	
 Function: Gets new genes from the gene_mapping table and inserts them in to the gene_events table  	
 Returns: Count of new genes  
 Args: Database handle object     
=cut

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

=head2 

 Title: get_changed_genes   	
 Usage: AnnotationEvents::get_changed_genes($dbh)   	
 Function: Gets exon_boundary exon_number CDS_change from the gene_mapping table and inserts them as change_gene in to the gene_events table  	
 Returns: Count of changed genes 
 Args: Database handle object	     
=cut

sub get_changed_genes {
	my($dbh) = @_;
	
	my $gene_count=0;
	my $changed_gene_insert_sth = $dbh->prepare(get_sql('gene_events'));
	
	my @maptypes = qw(exon_boundary exon_number CDS_change);
	my @changed_genes;
	
	foreach my $maptype (@maptypes){
		push @changed_genes, @{GeneMapping::get_gene_mappings_by_maptype($dbh,$maptype);};	
	}
	
	foreach my $row (@changed_genes){
			my $vb_gene_id  = $row->[1];
		    my $cap_gene_id = $row->[0];
			$changed_gene_insert_sth->execute($vb_gene_id,$cap_gene_id,'change_gene');
			$gene_count++;
	}
	
	return $gene_count;
}

=head2
 Title: get_lost_iso_form   	
 Usage: AnnotationEvents::get_splits($dbh)   	
 Function: Gets lost_iso_form from the gene_mapping table and inserts them in to the gene_events table 	
 Returns: Total count of genes with lost_iso_form
 Args: Database handle object 
=cut

sub get_lost_iso_form {
	my($dbh) = @_;
	
	my $gene_count=0;
	my $changed_gene_insert_sth = $dbh->prepare(get_sql('gene_events'));
	
	my @maptypes = qw(lost_iso_form);
	my @changed_genes;
	
	foreach my $maptype (@maptypes){
		push @changed_genes, @{GeneMapping::get_gene_mappings_by_maptype($dbh,$maptype);};	
	}
	
	foreach my $row (@changed_genes){
			my $vb_gene_id  = $row->[1];
		    my $cap_gene_id = $row->[0];
			$changed_gene_insert_sth->execute($vb_gene_id,$cap_gene_id,'lost_iso_form');
			$gene_count++;
	}
	
	return $gene_count;
}

=head2
 Title: get_gain_iso_form   	
 Usage: AnnotationEvents::get_gain_iso_form($dbh)   	
 Function: Gets gain_iso_form from the gene_mapping table and inserts them in to the gene_events table 	
 Returns: Total count of genes with gain_iso_form
 Args: Database handle object 
=cut

sub get_gain_iso_form {
	my($dbh) = @_;
	
	my $gene_count=0;
	my $changed_gene_insert_sth = $dbh->prepare(get_sql('gene_events'));
	
	my @maptypes = qw(gain_iso_form);
	my @changed_genes;
	
	foreach my $maptype (@maptypes){
		push @changed_genes, @{GeneMapping::get_gene_mappings_by_maptype($dbh,$maptype);};	
	}
	
	foreach my $row (@changed_genes){
			my $vb_gene_id  = $row->[1];
		    my $cap_gene_id = $row->[0];
			$changed_gene_insert_sth->execute($vb_gene_id,$cap_gene_id,'gain_iso_form');
			$gene_count++;
	}
	
	return $gene_count;
}

=head2 

 Title: get_splits   	
 Usage: AnnotationEvents::get_splits($dbh)   	
 Function: Gets split genes from the gene_mapping table and inserts them in to the gene_events table 	
 Returns: Total count of splits and count of core and cap genes involved
 Args: Database handle object    
=cut

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

=head2 

 Title: get_merge   	
 Usage: AnnotationEvents::get_merge($dbh)   	
 Function: Gets merged genes from the gene_mapping table and inserts them in to the gene_events table 	
 Returns: Total count of merges and count of core and cap genes involved
 Args: Database handle object     
=cut

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
