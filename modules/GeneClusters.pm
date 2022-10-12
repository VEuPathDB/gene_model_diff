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

GeneClusters

=head1 SYNOPSIS
	
	use GeneClusters;
	
	GeneClusters::work_out_gene_clusters($dbh);

=head1 DESCRIPTION

This module groups core and cap genes together based on exon overlap.
This module is the interface to the gene_clusters table.

=head1 Author

	Mikkel B Christensen

=head1 METHODS

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

=head2 work_out_gene_clusters

 Title: work_out_gene_clusters
 Usage: GeneClusters::work_out_gene_clusters($dbh)
 Function: Group genes based exon overlap
 Returns: Nothing
 Args: Database handle object
=cut

sub work_out_gene_clusters {
  my ($dbh) = @_;
  my %added_cap_genes;
  my $cap_gene_model_ref = GeneModel::get_distinct_id_by_source($dbh, 'gene', 'cap');
  foreach my $cap_gene_id_ref (@{$cap_gene_model_ref}) {
    if (!exists $added_cap_genes{$cap_gene_id_ref->[0]}) {
      my %gene_cluster;
      _recursive_cluster_genes($dbh, $cap_gene_id_ref, \%gene_cluster);
      _insert_gene_cluster($dbh, \%gene_cluster);
      foreach my $cap_gene_id (keys %{$gene_cluster{cap}}) {
        $added_cap_genes{$cap_gene_id} = 1;
      }
    }
  }
}

=head2 calculate_cluster_summary

 Title: calculate_cluster_summary
 Usage: GeneClusters::calculate_cluster_summary($dbh)
 Function: count the number of genes and transcript in each cluster
 Returns: Nothing
 Args: Database handle object
=cut

sub calculate_cluster_summary {
  my ($dbh) = @_;
  my $clusters_ref = get_gene_clusters($dbh);
  my %clusters;
  foreach my $cluster_row (@{$clusters_ref}) {
    my ($cluster_id, $gene_id, $source) = @{$cluster_row};

    my $transcripts =
      GeneModel::get_distinct_id_by_id($dbh, 'gene', $gene_id, 'transcript', $source, 1);

    $clusters{$cluster_id}{$source}{$gene_id} = scalar(@{$transcripts});
  }

  foreach my $cluster_id (keys %clusters) {
    my ($cap_gene_count, $cap_trans_count, $vb_gene_count, $vb_trans_count);
    foreach my $source (keys %{$clusters{$cluster_id}}) {
      my $gene_count;
      my $transcript_count;
      foreach my $gene_id (keys %{$clusters{$cluster_id}{$source}}) {
        $gene_count++;
        $transcript_count += $clusters{$cluster_id}{$source}{$gene_id};
      }
      if ($source eq 'cap') {
        $cap_gene_count  = $gene_count;
        $cap_trans_count = $transcript_count;
      } elsif ($source eq 'vb') {
        $vb_gene_count  = $gene_count;
        $vb_trans_count = $transcript_count;
      }
    }

    _insert_cluster_summary($dbh, $cluster_id, $cap_gene_count, $cap_trans_count, $vb_gene_count,
      $vb_trans_count);
  }
}

=head2 get_gene_clusters

 Title:	get_gene_clusters
 Usage:	GeneClusters::get_gene_clusters($dbh,$array_ref)
 Function: select all cluster data for error range. 
 Returns: Array ref of cluster data
 Args: Database handle object, Array ref of lower and upper error limit
=cut

sub get_gene_clusters {
  my ($dbh, $error_codes) = @_;
  my $errorLog = Log::Log4perl->get_logger("errorlogger");
  if (!defined $error_codes) {
    $errorLog->error(
      "GeneClusters::get_gene_clusters needs lower and upper limit error code. default -100:100");
    my @array = (-100, 100);
    $error_codes = \@array;
  } elsif (scalar @{$error_codes} != 2) {
    $errorLog->error(
      "GeneClusters::get_gene_clusters needs lower and upper limit error code. default -100:100");
    @{$error_codes} = ();
    @{$error_codes} = (-100, 100);
  }

  my $sql =
    "select gene_cluster_id,gene_id,source from gene_clusters where error_code between ? and ? ;";
  my $array_ref = _execute_sql($dbh, $sql, $error_codes);
  return $array_ref;
}

=head2 get_distinct_cluster_ids

 Title: get_distinct_cluster_ids
 Usage: GeneClusters::get_distinct_cluster_ids($dbh,$array_ref)
 Function: get id for all clusters in error range.
 Returns: Array_ref of cluster IDs
 Args: Database handle object, Array ref of lower and upper error limit
=cut

sub get_distinct_cluster_ids {
  my ($dbh, $error_codes) = @_;
  my $errorLog = Log::Log4perl->get_logger("errorlogger");
  if (!defined $error_codes) {
    $errorLog->error(
"GeneClusters::get_distinct_cluster_ids needs lower and upper limit error code. default -100:100"
    );
    my @array = (-100, 100);
    $error_codes = \@array;
  } elsif (scalar @{$error_codes} != 2) {
    $errorLog->error(
"GeneClusters::get_distinct_cluster_ids needs lower and upper limit error code. default -100:100"
    );
    @{$error_codes} = ();
    @{$error_codes} = (-100, 100);
  }

  my $sql =
"select distinct(gene_cluster_id) from gene_clusters where source = 'cap' group by gene_cluster_id having min(error_code) >= ? and max(error_code) <= ?;";
  my $array_ref = _execute_sql($dbh, $sql, $error_codes);

  my @ids;
  foreach my $row_ref (@{$array_ref}) {
    push @ids, $row_ref->[0];
  }

  return \@ids;
}

=head2 get_gene_cluster_by_id

 Title: get_gene_cluster_by_id
 Usage: GeneClusters::get_gene_cluster_by_id($dbh,$cluster_id)
 Function: get data for a single cluster
 Returns: Array_ref with cluster data
 Args: Database handle object, cluster ID 
=cut

sub get_gene_cluster_by_id {
  my ($dbh, $gene_cluster_id) = @_;
  my $sql = "select gene_id,source from gene_clusters where gene_cluster_id = $gene_cluster_id;";
  my $array_ref = _submit_sql($dbh, $sql);
  return $array_ref;
}

=head2 get_cluster_summary

 Title: get_cluster_summary
 Usage: GeneClusters::get_cluster_summary($dbh)
 Function: get summary for all clusters
 Returns: Array ref
 Args: Database handle object
=cut

sub get_cluster_summary {
  my ($dbh) = @_;
  my $sql = "select  gene_cluster_id,
			   cap_gene_count,
			   cap_transcript_count,
			   vb_gene_count,
			   vb_transcript_count
			   from cluster_summary;";
  my $array_ref = _submit_sql($dbh, $sql);
  return $array_ref;
}

=head2 get_cluster_summary_by_id

 Title: get_cluster_summary_by_id
 Usage: GeneClusters::get_cluster_summary_by_id()
 Function: get summary for a single cluster
 Returns: Array_ref, Cluster ID
 Args: Database handle object, cluster ID
=cut

sub get_cluster_summary_by_id {
  my ($dbh, $cluster_id) = @_;
  my $sql = "select cap_gene_count,
			   cap_transcript_count,
			   vb_gene_count,
			   vb_transcript_count,
			   cap_max_error,
			   vb_max_error
			   from cluster_summary
			   where gene_cluster_id = $cluster_id;";
  my $array_ref = _submit_sql($dbh, $sql);
  return $array_ref;
}

=head2 get_max_error

 Title: get_max_error
 Usage: GeneClusters::get_max_error()
 Function: get the max error code for cap and vb genes respectively
 Returns: Array_ref, Cluster ID
 Args: Database handle object, cluster ID
=cut

sub get_max_error {
  my ($dbh, $cluster_id) = @_;
  my $sql       = "select max(error_code) from gene_clusters where source = ?;";
  my @cap_param = ('cap');
  my @vb_param  = ('vb');

  my $cap_array = _execute_sql($dbh, $sql, \@cap_param);
  my $vb_array  = _execute_sql($dbh, $sql, \@vb_param);

  my $max_cap = $cap_array->[0]->[0];
  my $max_vb  = $vb_array->[0]->[0];

  return ($max_cap, $max_vb);
}

sub _recursive_cluster_genes {
  my ($dbh, $cap_gene_array, $gene_cluster) = @_;
  if (scalar @{$cap_gene_array} == 0) {
    return 1;
  }

  my $cap_gene_id = shift @{$cap_gene_array};
  if (!exists $gene_cluster->{cap}{$cap_gene_id}) {
    $gene_cluster->{cap}{$cap_gene_id} = 1;
    my $vb_gene_ids = _get_gene_via_mapped_exons($dbh, $cap_gene_id, 'cap', 'vb');

    foreach my $vb_gene_id (@{$vb_gene_ids}) {
      if (!exists $gene_cluster->{vb}{$vb_gene_id}) {
        $gene_cluster->{vb}{$vb_gene_id} = 1;
        my $cap_gene_ids = _get_gene_via_mapped_exons($dbh, $vb_gene_id, 'vb', 'cap');
        push(@{$cap_gene_array}, @{$cap_gene_ids});
      }
    }
  }
  _recursive_cluster_genes($dbh, $cap_gene_array, $gene_cluster);
}

sub _get_gene_via_mapped_exons {
  my ($dbh, $query_gene_id, $query_source, $return_source) = @_;
  my %gene_id_duplicate;
  my @return_gene_ids;
  my @query_exon_ids =
    @{GeneModel::get_gene_model_by_id($dbh, 'gene', $query_gene_id, $query_source, 'exon', 1)};
  foreach my $exon_id (@query_exon_ids) {
    my $return_exon_ids = ExonMapping::get_exon_mappings($dbh, $query_source, $exon_id);
    foreach my $row (@{$return_exon_ids}) {
      my $return_exon_id = @{$row}[1];
      next unless defined($return_exon_id);
      my $return_gene_id =
        GeneModel::get_distinct_id_by_id($dbh, 'exon', $return_exon_id, 'gene', $return_source);
      if (defined($return_gene_id) and !exists $gene_id_duplicate{$return_gene_id}) {
        $gene_id_duplicate{$return_gene_id} = 1;
        push @return_gene_ids, $return_gene_id;
      }
    }
  }

  @return_gene_ids =
    _filter_overlapping_cds($dbh, $query_gene_id, \@return_gene_ids, $query_source, $return_source);

  return \@return_gene_ids;
}

# Filter out all gene pairs that do not have any overlapping CDSs
sub _filter_overlapping_cds {
  my ($dbh, $query_gene_id, $return_gene_ids, $query_source, $return_source) = @_;

  # Get all CDSs for the query gene
  my $query_transcripts =
    GeneModel::get_transcripts_by_gene_id($dbh, $query_gene_id, $query_source);
  my @query_cdss = ();
  for my $query_tr (@$query_transcripts) {
    my $qcdss = CDS::get_all_cds_pos_by_parent_id($dbh, $query_tr);
    push @query_cdss, @$qcdss;
  }
  @query_cdss = sort { $a->[0] <=> $b->[0] or $a->[1] <=> $b->[1] } @query_cdss;

  my @ok_return_gene_ids;
  for my $return_gene_id (@$return_gene_ids) {

    # Get all CDSs for that return gene
    my $return_transcripts =
      GeneModel::get_transcripts_by_gene_id($dbh, $return_gene_id, $return_source);

    my @return_cdss;
    for my $transcript (@$return_transcripts) {
      my $rcdss = CDS::get_all_cds_pos_by_parent_id($dbh, $transcript);
      push @return_cdss, @$rcdss;
    }
    @return_cdss = sort { $a->[0] <=> $b->[0] or $a->[1] <=> $b->[1] } @return_cdss;

    # See if the CDSs overlap
    # Only keep the genes that overlap even a little with a CDS
    if (_overlapping_regions(\@query_cdss, \@return_cdss)) {
      push @ok_return_gene_ids, $return_gene_id;
    }
  }

  return @ok_return_gene_ids;
}

# Very basic comparison of 2 sets of ordered arrays of regions
# Each an array of 2 values: start and end (integers)
sub _overlapping_regions {
  my ($regions1, $regions2) = @_;

REG1: for my $reg1 (@$regions1) {
    my ($start1, $end1) = @$reg1;

  REG2: for my $reg2 (@$regions2) {
      my ($start2, $end2) = @$reg2;

      next REG1 if $start2 > $end1;
      next REG2 if $start1 > $end2;

      # Overlap!
      return 1;
    }
  }

  return 0;
}

sub _insert_gene_cluster {
  my ($dbh, $gene_cluster) = @_;
  my $insert_sql = "insert gene_clusters(
			   gene_cluster_id, 
			   gene_id,
			   source,
			   error_code)
			   select ?,?,?,?;";
  my $insert_sth = $dbh->prepare($insert_sql);

  #get cluster ID!
  my $get_last_cluster_id_sql = "select max(gene_cluster_id) from gene_clusters;";
  my $last_cluster_ref        = _submit_sql($dbh, $get_last_cluster_id_sql);
  my $last_cluster_id = defined($last_cluster_ref->[0]->[0]) ? $last_cluster_ref->[0]->[0] : 0;
  my $gene_cluster_id = $last_cluster_id + 1;
  foreach my $gene_id (keys %{$gene_cluster->{cap}}) {
    my $error_code = GeneModel::get_error_code_by_gene_id($dbh, $gene_id, 'cap');
    $insert_sth->execute($gene_cluster_id, $gene_id, 'cap', $error_code);
  }

  foreach my $gene_id (keys %{$gene_cluster->{vb}}) {
    my $error_code = GeneModel::get_error_code_by_gene_id($dbh, $gene_id, 'vb');
    $insert_sth->execute($gene_cluster_id, $gene_id, 'vb', $error_code);
  }

}

sub _insert_cluster_summary {
  my ($dbh, $gene_cluster_id, $cap_gene_count, $cap_trans_count, $vb_gene_count, $vb_trans_count) =
    @_;

  if (!defined($vb_gene_count)) {
    $vb_gene_count  = 0;
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
  $insert_sth->execute($gene_cluster_id, $cap_gene_count, $cap_trans_count, $vb_gene_count,
    $vb_trans_count);
}

sub _execute_sql {
  my ($dbh, $sql, $params) = @_;
  my $sth = $dbh->prepare($sql);
  $sth->execute(@{$params});
  return $sth->fetchall_arrayref();
}

sub _submit_sql {
  my ($dbh, $sql) = @_;
  my $array_ref = $dbh->selectall_arrayref($sql);
  return $array_ref;
}

1;

