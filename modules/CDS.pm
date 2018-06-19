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

=head1 CONTACT
	
	Please email comments or questions to info@vectorbase.org
	
=cut

=head1 NAME

CDS

=head1 SYNOPSIS
	
	use CDS;
	
	my($start,$end) = CDS::get_cds_pos_by_parent_id($database_handle,$parent_id); 

=head1 DESCRIPTION

This module is the interface to the CDS table which contains genomic coordinates and MD5 checksum for all CDS.
The module contains the SQL to extract start and end position and MD5 checksum for each CDS using the ID for the parent mRNA.

=head1 Author

	Mikkel B Christensen

=head1 METHODS

=cut

package CDS;

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use DBI;

=head2 get_all_cds_by_parent_id

 Title:    get_all_cds_by_parent_id	
 Usage:    CDS::get_cds_pos_by_parent_id($dbh,$parent_id);
 Function: get start and end position for CDS. 	
 Returns:  list (start,end). 	
 Args:     Database handle object,parent mRNA ID. 
=cut 

sub get_cds_pos_by_parent_id{
	my($dbh,$parent_id) = @_;
	
	my $sql = "select start,end from cds where cds_parent_id = \'$parent_id\';";
	my $array_ref = $dbh->selectall_arrayref($sql);
	my $start = $array_ref->[0]->[0];
	my $end   = $array_ref->[0]->[1];
	return($start,$end);	
};

=head2 get_cds_checksum_by_parent_id

 Title:    get_cds_checksum_by_parent_id	
 Usage:    CDS::get_cds_checksum_parent_id($dbh,$parent_id);
 Function: get checksum for CDS. 	
 Returns:  string.
 Args:     Database handle object,parent mRNA ID. 
=cut 

sub get_cds_checksum_by_parent_id{
	my($dbh,$parent_id) = @_;
	
	my $sql = "select md5_checksum from cds where cds_parent_id = \'$parent_id\';";
	my $array_ref = $dbh->selectall_arrayref($sql);
	
	return $array_ref->[0]->[0];
}

=head2 get_cds_error_code_by_parent_id

 Title:    get_cds_error_code_by_parent_id	
 Usage:    CDS::get_cds_error_code_by_parent_id($dbh,$parent_id);
 Function: get cds_error_code for CDS. 	
 Returns:  string.
 Args:     Database handle object,parent mRNA ID. 
=cut 

sub get_cds_error_code_by_parent_id{
	my($dbh,$parent_id) = @_;
	
	my $sql = "select cds_error_code from cds where cds_parent_id = \'$parent_id\';";
	my $array_ref = $dbh->selectall_arrayref($sql);
	
	return $array_ref->[0]->[0];
}

1;