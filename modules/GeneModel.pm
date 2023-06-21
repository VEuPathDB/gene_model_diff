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

GeneModel

=head1 SYNOPSIS

	use GeneModel;
	
	GeneModel::get_gene_model_by_id($dbh,$type,$id,$source,$return_type,$force_array);

=head1 DESCRIPTION

This module is the interface to the gene_model table.

=head1 Author

	Mikkel B Christensen

=head1 METHODS

=cut

package GeneModel;
use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use DBI;

# Define supported biotypes
require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT_OK = qw(%BIOTYPE);
my %known_biotypes = (
  gene => [qw(
    gene
    protein_coding_gene
    pseudogene
    ncRNA_gene
  )],
  transcript => [qw(
    mRNA
    transcript
    pseudogenic_transcript
    ncRNA
    lnc_RNA
    scRNA
    snRNA
    snoRNA
    rRNA
    tRNA
  )],
  sub_feature => [qw(
    exon
    CDS
    non_canonical_five_prime_splice_site
    non_canonical_three_prime_splice_site
  )],
  skip => [qw(
    region
    five_prime_UTR
    three_prime_UTR

  )]
);

our %BIOTYPE = ();
for my $name (keys %known_biotypes) {
  $BIOTYPE{$name} = {map { $_ => 1 } @{$known_biotypes{$name}}};
}

my %types = (
  "gene"       => "gene_id",
  "transcript" => "transcript_id",
  "exon"       => "exon_id"
);

my %return_types = (
  "gene"       => "2",
  "transcript" => "1",
  "exon"       => "0",
  "gene_model" => "0..2"
);

=head2 get_transcripts_by_gene_id

 Title: get_transcripts_by_gene_id
 Usage: GeneModel::get_transcripts_by_gene_id($dbh,$gene_id,$source)
 Function: Get Transcript IDs for a gene ID
 Returns: Array ref of transcript IDs
 Args: Database handle object, Gene ID, Source of gene 
=cut

sub get_transcripts_by_gene_id {
  my ($dbh, $gene_id, $source) = @_;
  my $array_ref =
    GeneModel::get_distinct_id_by_id($dbh, 'gene', $gene_id, 'transcript', $source, 1);
  my @transcript_ids;
  foreach my $row (@{$array_ref}) {
    push @transcript_ids, $row->[0];
  }
  return \@transcript_ids;

}

=head2 get_exons_by_transcript_id

 Title: get_exons_by_transcript_id
 Usage: GeneModel::get_exons_by_transcript_id($dbh,$transcript_id,$source)
 Function: Get exon IDs for transcript ID
 Returns: Array ref exon IDs
 Args: Database handle object, transcript ID, Source of transcript.
=cut

sub get_exons_by_transcript_id {
  my ($dbh, $transcript_id, $source) = @_;
  my $gene_model_ref = GeneModel::get_gene_model_by_id($dbh, 'transcript', $transcript_id, $source);
  my @exons;
  foreach my $row (@{$gene_model_ref}) {
    push @exons, $row->[0];
  }
  return \@exons;
}

=head2 

 Title: get_gene_model_by_id
 Usage: GeneModel::get_gene_model_by_id($dbh,$type,$id,$source,$return_type,$force_array)
 Function: Get gene_model by either Gene/Transcript/Exon ID
 Returns: Array_ref or single feature ID  
 Args: Database handle object, type of feature ID to search with,source,feature ID, feature ID to return, boolean to force array.  
=cut

sub get_gene_model_by_id {
  my ($dbh, $type, $id, $source, $return_type, $force_array) = @_;
  confess("missing parameter to get_gene_model_by_id") unless ($dbh and $type and $id);
  my $sql;
  if ($source) {
    $sql = "select exon_id,transcript_id,gene_id,source
			   	   from gene_model where $types{$type} = \'$id\' and source = \'$source\';";
  } else {
    $sql = "select exon_id,transcript_id,gene_id,source
			   	   from gene_model where $types{$type} = \'$id\';";
  }

  my $array_ref = _submit_sql($dbh, $sql);

  if (scalar @{$array_ref} == 1 and $return_type and $force_array) {
    my @IDs;
    foreach my $row (@{$array_ref}) {
      my $id = @{$row}[$return_types{$return_type}];
      push @IDs, $id;
    }
    return \@IDs;
  } elsif (scalar @{$array_ref} == 1 and $return_type) {
    return $array_ref->[0][$return_types{$return_type}];

  } elsif (scalar @{$array_ref} > 1 and $return_type) {
    my @IDs;
    foreach my $row (@{$array_ref}) {
      my $id = @{$row}[$return_types{$return_type}];
      push @IDs, $id;
    }
    return \@IDs;
  } else {
    return $array_ref;
  }
}

=head2 get_distinct_id_by_id

 Title: get_distinct_id_by_id
 Usage: GeneModel::get_distinct_id_by_id($dbh,$type,$id,$return_id,$source,$force_array)
 Function: Get Gene/Transcript/exon/ ID seaech whith Gene/Transcript/exon/ ID  
 Returns: Single feature ID
 Args: Database handle object, type of feature ID to search with,feature ID,feature ID to return,source,boolean to force array. 
=cut

sub get_distinct_id_by_id {
  my ($dbh, $type, $id, $return_id, $source, $force_array) = @_;
  confess("missing parameter to get_distinct_id_by_id")
    unless ($dbh and $type and $id and $return_id and $source);
  my $sql =
"select distinct $types{$return_id} from gene_model where $types{$type} = \'$id\' and source = \'$source\';";
  my $array_ref = _submit_sql($dbh, $sql);

  if ($force_array) {
    return $array_ref;
  } else {
    return $array_ref->[0]->[0];
  }

}

=head2 get_distinct_id_by_source

 Title: get_distinct_id_by_source
 Usage: GeneModel::get_distinct_id_by_source($dbh,$type,$source)
 Function: get all Get Gene/Transcript/exon/ ID for a source
 Returns array ref of feature ID
 Args: Database handle object, feature type, source.
=cut

sub get_distinct_id_by_source {
  my ($dbh, $type, $source) = @_;

  my $sql       = "select distinct $types{$type} from gene_model where source = \'$source\';";
  my $array_ref = _submit_sql($dbh, $sql);
  return $array_ref;
}

=head2 get_error_code_by_gene_id

 Title: get_error_code_by_gene_id
 Usage: GeneModel::get_error_code_by_gene_id($dbh,$gene_id,$source)
 Function: get error code for gene_model
 Returns array_ref
 Args: Database handle object, Gene ID, Source
=cut

sub get_error_code_by_gene_id {
  my ($dbh, $gene_id, $source) = @_;
  my $sql =
"select distinct error_code from gene_model where gene_id = \'$gene_id\' and source = \'$source\';";
  my $array_ref  = _submit_sql($dbh, $sql);
  my $error_code = $array_ref->[0]->[0];
  return $error_code;
}

sub _submit_sql {
  my ($dbh, $sql) = @_;
  my $array_ref = $dbh->selectall_arrayref($sql);
  return $array_ref;
}

1;
