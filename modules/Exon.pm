
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

Exon

=head1 SYNOPSIS

  use Exon;
  
  my exon_overlap_Array_ref = Exon::get_overlapping_exons($dbh,'partial_3');
  
=head1 DESCRIPTION
  
This module is the interface to the exon table which contains genomic coordinates for all Exons.
The module contains the SQL to extract exons that overlap between core and cap geneset.

=head1 Author

  Mikkel B Christensen

=head1 METHODS

=cut

package Exon;

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;

=head2 get_overlapping_exons

 Title:    get_overlapping_exons  
 Usage:    Exon::get_overlapping_exons($dbh,'partial_3').  
 Function: Execute SQL on the database.   
 Returns:  Array refference of pairs of overlapping exon ids.   
 Args:     Database handle object,type of overlap, where type can be one of the following strings (identical|partial_3|partial_5|spanning|included|no_match_cap). 
=cut 

sub get_overlapping_exons {
  my ($dbh, $overlap_type) = @_;
  my $sql = _get_sql($overlap_type);

  #log error if no SQL
  unless ($sql) {
    confess("There is no SQL of type $overlap_type, therefore no Exons will be returned");
  }
  my $sth = $dbh->prepare($sql);

  $sth->execute();
  my $array_ref = $sth->fetchall_arrayref;

  return $array_ref;
}

sub _get_sql {
  my ($sql_name, $same_biotype) = @_;

  $same_biotype //= 1;
  my $same_biotype_sql = " and gene_cap.biotype = gene_vb.biotype";
  
  my $base_sql = "SELECT cap.exon_id, vb.exon_id
    FROM
    exon cap
      LEFT JOIN gene_model gene_cap ON (cap.exon_id = gene_cap.exon_id),
    exon vb
      LEFT JOIN gene_model gene_vb ON (vb.exon_id = gene_vb.exon_id)
    WHERE cap.source   = 'cap'
      and vb.source    = 'vb'
      and cap.scaffold = vb.scaffold
      and cap.strand   = vb.strand
  ";

  if ($sql_name eq 'identical') {
    my $sql = $base_sql . "
      and cap.start = vb.start
      and cap.end   = vb.end";
    $sql .= $same_biotype_sql if $same_biotype;
    return $sql;
  } elsif ($sql_name eq 'partial_3') {
    my $sql = $base_sql . "
      and cap.end  >= vb.start
      and cap.end  <= vb.end";
    $sql .= $same_biotype_sql if $same_biotype;
    return $sql;
  } elsif ($sql_name eq 'partial_5') {
    my $sql = $base_sql . "
      and cap.start  >= vb.start
      and cap.start  <= vb.end";
    $sql .= $same_biotype_sql if $same_biotype;
    return $sql;

  } elsif ($sql_name eq 'spanning') {
    my $sql = $base_sql . "
      and cap.start  <= vb.start
      and cap.end  >= vb.end";
    $sql .= $same_biotype_sql if $same_biotype;
    return $sql;

  } elsif ($sql_name eq 'included') {
    my $sql = $base_sql . "
      and cap.start  >= vb.start
      and cap.end  <= vb.end";
    $sql .= $same_biotype_sql if $same_biotype;
    return $sql;

  } elsif ($sql_name eq 'no_match_cap') {
    my $sql = "select exon_id, source from exon
           where source = 'cap'
           and  not exists (select * from exon_mappings where exon_id = cap_exon_id)";
    return $sql;
  }
}
1;
