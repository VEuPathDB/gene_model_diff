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

TranscriptMapping

=head1 SYNOPSIS

TranscriptMapping::resolve_transcript_mappings($dbh)

=head1 DESCRIPTION

In each gene cluster select the set of transcript pairs with least differences lowest rank.

=head1 Author

	Mikkel B Christensen

=head1 METHODS

=cut

package TranscriptMapping;
use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use TranscriptLinks;
use Data::Dumper;

=head2 resolve_transcript_mappings

 Title: resolve_transcript_mappings
 Usage: TranscriptMapping::resolve_transcript_mappings($dbh)
 Function: In each gene cluster select the set of transcript pairs with lowest rank
 Returns: Populates the the gene_mappings table
 Args: Database handle
=cut 

sub resolve_transcript_mappings {
  my ($dbh) = @_;
  my $select_transcript_by_rank_sql =
"select id,cap_transcript_id,vb_transcript_id,group_count from transcript_links where gene_cluster_id = ? and link_status = 'not_mapped' group by cap_transcript_id,vb_transcript_id having max(link_rank) = ?;";
  my $select_transcript_by_rank_sth = $dbh->prepare($select_transcript_by_rank_sql);
  my @error_limits                  = (0, 9);
  my $gene_cluster_ids              = GeneClusters::get_distinct_cluster_ids($dbh, \@error_limits);
  for my $cluster_id (@{$gene_cluster_ids}) {
    _loop_over_rank($dbh, $cluster_id, $select_transcript_by_rank_sth);
  }
  _copy_from_link_to_mapping_table($dbh);
}

=head2 get_all_transcript_mappings

 Title: get_all_transcript_mappings
 Usage: TranscriptMapping::get_all_transcript_mappings($dbh)
 Function: selects all transcript mappings
 Returns: Array ref
 Args: Database handle
=cut 

sub get_all_transcript_mappings {
  my ($dbh)     = @_;
  my $sql       = "select cap_trans_id,vb_trans_id,map_type from transcript_mappings;";
  my $array_ref = $dbh->selectall_arrayref($sql);
  return $array_ref;
}

=head2 get_all_transcript_mappings_by_id

 Title: get_all_transcript_mappings_by_id
 Usage: TranscriptMapping::get_all_transcript_mappings_by_id($dbh,$cluster_id)
 Function: selects transcripts mappings for a single cluster
 Returns: Array ref
 Args: Database handle
=cut 

sub get_all_transcript_mappings_by_id {
  my ($dbh, $cluster_id) = @_;
  my %attr;
  my @values;
  push @values, $cluster_id;
  my $sql =
"select cap_trans_id,vb_trans_id,map_type from transcript_mappings where gene_cluster_id = \'$cluster_id\';";
  my $array_ref = $dbh->selectall_arrayref($sql);
  return $array_ref;
}

sub _loop_over_rank {
  my ($dbh, $cluster_id, $select_transcript_by_rank_sth) = @_;
  for (my $rank = 1 ; $rank <= 5 ; $rank++) {
    $select_transcript_by_rank_sth->execute($cluster_id, $rank);
    my $array_ref = $select_transcript_by_rank_sth->fetchall_arrayref;
    next unless (scalar @{$array_ref} > 0);
    my @transcript_pairs_to_keep;
    _get_transcript_pair_with_min_count($array_ref, \@transcript_pairs_to_keep);
    _update_transcript_link_table($dbh, \@transcript_pairs_to_keep);
  }

  if (_find_unmapped_transcript($dbh, $cluster_id)) {
    _loop_over_rank($dbh, $cluster_id, $select_transcript_by_rank_sth);
  }
}

sub _get_transcript_pair_with_min_count {
  my ($array_ref, $transcript_pairs_to_keep) = @_;
  if (scalar @{$array_ref} == 0) {
    return;
  } elsif (scalar @{$array_ref} == 1) {
    push @{$transcript_pairs_to_keep}, @{$array_ref};
    return;
  }

  my $uniq_pairs = _get_lowest_count($array_ref);
  push @{$transcript_pairs_to_keep}, @{$uniq_pairs};
  _remove_redundant_transcript($array_ref, $uniq_pairs);
  _get_transcript_pair_with_min_count($array_ref, $transcript_pairs_to_keep);
}

sub _get_lowest_count {
  my ($array_ref) = @_;

  my @sorted_array = sort { $a->[3] <=> $b->[3] } @{$array_ref};
  my $row          = shift @sorted_array;
  my $lowest_count = $row->[3];

  push(my (@least_count), $row);
  foreach my $row (@sorted_array) {
    my $count = $row->[3];
    if ($count == $lowest_count) {
      push(my (@least_count), $row);
    } else {
      last;
    }
  }

  if (scalar @least_count == 1) {
    return \@least_count;
  } else {
    my $uniq_pairs = _is_transcript_uniq(\@least_count);
    return $uniq_pairs;
  }
}

sub _is_transcript_uniq {
  my ($least_count) = @_;

  my %cap_ids;
  my %vb_ids;
  my @uniq_pair;
  foreach my $transcript_pair (@{$least_count}) {
    my ($id, $cap_transcript_id, $vb_transcript_id, $group_count) = @{$transcript_pair};

    if ((exists $cap_ids{$cap_transcript_id}) or (exists $vb_ids{$vb_transcript_id})) {

      #write to log!
    } else {
      $cap_ids{$cap_transcript_id} = 1;
      $vb_ids{$vb_transcript_id}   = 1;
      push @uniq_pair, $transcript_pair;
    }
  }
  return \@uniq_pair;
}

sub _remove_redundant_transcript {
  my ($array_ref, $uniq_pairs) = @_;

  my %cap_ids;
  my %vb_ids;

  foreach my $uniq_pair (@{$uniq_pairs}) {
    my ($id, $cap_transcript_id, $vb_transcript_id, $group_count) = @{$uniq_pair};
    $cap_ids{$cap_transcript_id} = 1;
    $vb_ids{$vb_transcript_id}   = 1;
  }

  my $index = 0;
  foreach my $transcript_pair (@{$array_ref}) {
    my ($id, $cap_transcript_id, $vb_transcript_id, $group_count) = @{$transcript_pair};
    if ((exists $cap_ids{$cap_transcript_id}) or (exists $vb_ids{$vb_transcript_id})) {
      splice(@{$array_ref}, $index, 1);
    } else {
      $index++;
    }
  }
}

sub _find_unmapped_transcript {
  my ($dbh, $cluster_id) = @_;

  my $cap_un_mapped_sql = "select distinct t1.cap_transcript_id 
	                     from transcript_links t1 
	                     where t1.gene_cluster_id = ? 
	                     and t1.link_status = 'unavailable' 
	                     and not exists (select * from transcript_links t2 where t1.cap_transcript_id = t2.cap_transcript_id and t2.gene_cluster_id = ? and t2.link_status = 'mapped');";

  my $vb_un_mapped_sql = "select distinct t1.vb_transcript_id 
	                     from transcript_links t1 
	                     where t1.gene_cluster_id = ? 
	                     and t1.link_status = 'unavailable' 
	                     and not exists (select * from transcript_links t2 where t1.vb_transcript_id = t2.vb_transcript_id and t2.gene_cluster_id = ? and t2.link_status = 'mapped');";

  my $cap_un_mapped_sth = $dbh->prepare($cap_un_mapped_sql);
  my $vb_un_mapped_sth  = $dbh->prepare($vb_un_mapped_sql);
  $cap_un_mapped_sth->execute($cluster_id, $cluster_id);
  $vb_un_mapped_sth->execute($cluster_id, $cluster_id);
  my $cap_array_ref = $cap_un_mapped_sth->fetchall_arrayref;
  my $vb_array_ref  = $vb_un_mapped_sth->fetchall_arrayref;

  if (scalar @{$cap_array_ref}) {
    _reset_unmapped_transcripts($dbh, $cap_array_ref, 'cap');
    return 1;
  } elsif (scalar @{$vb_array_ref}) {
    _reset_unmapped_transcripts($dbh, $vb_array_ref, 'vb');
    return 1;
  } else {
    return 0;
  }

}

sub _reset_unmapped_transcripts {
  my ($dbh, $array_ref, $source) = @_;

  my $cap_sql =
    "update transcript_links set link_status = 'not_mapped' where cap_transcript_id = ?";
  my $vb_sql  = "update transcript_links set link_status = 'not_mapped' where vb_transcript_id = ?";
  my $cap_sth = $dbh->prepare($cap_sql);
  my $vb_sth  = $dbh->prepare($vb_sql);

  foreach my $row (@{$array_ref}) {
    my $id = $row->[0];
    if ($source eq 'cap') {
      $cap_sth->execute($id);
    } elsif ($source eq 'vb') {
      $vb_sth->execute($id);
    }
  }
}

sub _update_transcript_link_table {
  my ($dbh, $array_ref) = @_;
  foreach my $uniq_transcript_pair (@{$array_ref}) {
    my ($id, $cap_transcript_id, $vb_transcript_id, $count) = @{$uniq_transcript_pair};
    my $set_row_mapped_sql =
      "update transcript_links set link_status = 'mapped' where id = \'$id\';";
    $dbh->do($set_row_mapped_sql);
    my $set_row_unavailable_sql =
"update transcript_links set link_status = 'unavailable' where link_status = 'not_mapped' and (cap_transcript_id = \'$cap_transcript_id\' or vb_transcript_id = \'$vb_transcript_id\');";
    $dbh->do($set_row_unavailable_sql);
  }
}

sub _copy_from_link_to_mapping_table {
  my ($dbh) = @_;

  my $link_table_sql =
"select gene_cluster_id, cap_transcript_id,vb_transcript_id,link_group from transcript_links where link_status = 'mapped';";
  my $insert_into_transcript_mapping =
    "insert transcript_mappings(gene_cluster_id,cap_trans_id,vb_trans_id,map_type) select ?,?,?,?;";
  my $insert_into_transcript_mapping_sth = $dbh->prepare($insert_into_transcript_mapping);
  my $array_ref                          = $dbh->selectall_arrayref($link_table_sql);

  foreach my $row (@{$array_ref}) {
    my ($gene_cluster_id, $cap_transcript_id, $vb_transcript_id, $map_type) = @{$row};

    $insert_into_transcript_mapping_sth->execute($gene_cluster_id, $cap_transcript_id,
      $vb_transcript_id, $map_type);
  }
}

1;
