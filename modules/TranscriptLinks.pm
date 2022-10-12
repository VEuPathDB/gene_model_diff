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

TranscriptLinks

=head1 SYNOPSIS

TranscriptLinks::work_out_transcript_links($dbh).

=head1 DESCRIPTION

To determine if a community annotated transcript (cap) has changed it must be compared to its parallel core transcript.
This is done by making all Cap core transcript pairs and assign a score based on exon and CDS similarity.
The transcript pairs are ranked and inserted into the transcript_links table. 

=head1 Author

	Mikkel B Christensen

=head1 METHODS

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
use CDS;

=head2 work_out_transcript_links

 Title:    work_out_transcript_links	
 Usage:    TranscriptLinks::work_out_transcript_links().	
 Function: Foreach cluster of genes, pair the transcript from the two genesets
 Returns:  Populates the transcript_links table. 	
 Args:     Database handle
=cut 

sub work_out_transcript_links {
  my ($dbh)                       = @_;
  my $errorLog                    = Log::Log4perl->get_logger("errorlogger");
  my $insert_transcript_links_sth = $dbh->prepare(_get_sql('insert_transcript_links'));
  my @error_limits                = (0, 9);
  my $gene_cluster_ids            = GeneClusters::get_distinct_cluster_ids($dbh, \@error_limits);
  foreach my $cluster_id (@{$gene_cluster_ids}) {
    my $gene_cluster     = GeneClusters::get_gene_cluster_by_id($dbh, $cluster_id);
    my $transcript_pairs = _make_transcript_pairs($dbh, $gene_cluster);
    foreach my $transcript_pair (@{$transcript_pairs}) {
      my $ranks = _rank_transcript_pair($dbh, $transcript_pair);
      my ($cap_id, $vb_id) = split /:/, $transcript_pair;
      unless ($ranks) {
        $errorLog->error("No rank was returned for transcript pair: $transcript_pair");
        next;
      }
      foreach my $group (keys %{$ranks}) {
        $insert_transcript_links_sth->execute(
          $cluster_id, $cap_id, $vb_id, $group,
          $ranks->{$group}->{rank},
          $ranks->{$group}->{count}, 'not_mapped'
        );
      }
    }
  }
}

sub _make_transcript_pairs {
  my ($dbh, $gene_cluster) = @_;
  my @cap;
  my @vb;
  my @pairs;
  foreach my $gene (@{$gene_cluster}) {
    my ($gene_id, $source) = @{$gene};
    my $transcripts = GeneModel::get_transcripts_by_gene_id($dbh, $gene_id, $source);
    if ($source eq 'cap') {
      push @cap, @{$transcripts};
    } elsif ($source eq 'vb') {
      push @vb, @{$transcripts};
    } else {
      croak("unknown source $source");
    }
  }

  foreach my $cap_id (@cap) {
    foreach my $vb_id (@vb) {
      push @pairs, "$cap_id:$vb_id";
    }
  }
  return \@pairs;

}

sub _rank_transcript_pair {
  my ($dbh, $transcript_pair) = @_;
  my %ranks = (
    identical     => 1,
    exon_boundary => 2,
    exon_number   => 3,
    CDS_change    => 4,
    CDS_error     => 5
  );

  my %groups = (
    identical   => 'identical',
    included    => 'exon_boundary',
    partial_5   => 'exon_boundary',
    partial_3   => 'exon_boundary',
    spanning    => 'exon_boundary',
    add_exon    => 'exon_number',
    delete_exon => 'exon_number',
    CDS_change  => 'CDS_change',
    CDS_error   => 'CDS_error'
  );
  my %maptype_count;
  my ($cap_transcript_id, $vb_transcript_id) = split /:/, $transcript_pair;
  my $cap_exons = GeneModel::get_exons_by_transcript_id($dbh, $cap_transcript_id, 'cap');
  my $vb_exons  = GeneModel::get_exons_by_transcript_id($dbh, $vb_transcript_id,  'vb');
  my ($exon_maptype_array, $no_exon_maptype_array) = _compare_exons($dbh, $cap_exons, $vb_exons);

  if (scalar @{$exon_maptype_array}) {
    my $cds_status = _compare_cds($dbh, $cap_transcript_id, $vb_transcript_id);

    foreach my $map_type (@{$exon_maptype_array}) {
      my $group = $groups{$map_type};
      $maptype_count{$group}{count}++;
    }

    foreach my $map_type (@{$no_exon_maptype_array}) {
      my $group = $groups{$map_type};
      $maptype_count{$group}{count}++;
    }

    if ($cds_status) {
      $maptype_count{$cds_status}{count} = 1;
    }

    foreach my $group (keys %maptype_count) {
      $maptype_count{$group}{rank} = $ranks{$group};
    }
    return \%maptype_count;
  } else {
    return undef;
  }

}

sub _compare_exons {
  my ($dbh, $cap_exons, $vb_exons) = @_;
  my %cap_exon_ids = map { $_ => 1 } @{$cap_exons};
  my %vb_exon_ids  = map { $_ => 1 } @{$vb_exons};
  my @exon_overlap_map_types;
  my @no_overlap_map_types;
  foreach my $cap_exon (keys %cap_exon_ids) {
    foreach my $vb_exon (keys %vb_exon_ids) {
      my $map_type = ExonMapping::get_map_type_for_exon_pair($dbh, $cap_exon, $vb_exon);
      if ($map_type) {
        push @exon_overlap_map_types, $map_type;
        delete $cap_exon_ids{$cap_exon};
        delete $vb_exon_ids{$vb_exon};
      }
    }
  }

  foreach my $cap_exon (keys %cap_exon_ids) {
    push @no_overlap_map_types, 'add_exon';
  }
  foreach my $vb_exon (keys %vb_exon_ids) {
    push @no_overlap_map_types, 'delete_exon';
  }

  return \@exon_overlap_map_types, \@no_overlap_map_types;
}

sub _compare_cds {
  my ($dbh, $cap_transcript_id, $vb_transcript_id) = @_;
  my ($cap_cds_checksum) = CDS::get_cds_checksum_by_parent_id($dbh, $cap_transcript_id);
  my ($cap_error_code)   = CDS::get_cds_error_code_by_parent_id($dbh, $cap_transcript_id);
  my ($vb_cds_checksum)  = CDS::get_cds_checksum_by_parent_id($dbh, $vb_transcript_id);
  my ($vb_error_code)    = CDS::get_cds_error_code_by_parent_id($dbh, $vb_transcript_id);

  if ($cap_error_code > $vb_error_code) {
    return 'CDS_error';
  }

  if ($cap_cds_checksum ne $vb_cds_checksum) {
    return 'CDS_change';
  } else {
    return 0;
  }

}

sub _get_sql {

  my ($sql_name) = @_;

  if ($sql_name eq 'insert_transcript_links') {
    my $insert_transcript_links =
"insert transcript_links(gene_cluster_id,cap_transcript_id,vb_transcript_id,link_group,link_rank,group_count,link_status) select ?,?,?,?,?,?,?";
    return $insert_transcript_links;
  }
}
1;
