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

package CDS;

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use DBI;

sub get_all_cds_by_parent_id{
	my($dbh,$parent_id) = @_;
	
	my $sql = "select start,end from cds where cds_parent_id = \'$parent_id\';";
	my $array_ref = $dbh->selectall_arrayref($sql);
	return $array_ref;	
};

sub get_cds_start_by_parent_id{
	my($dbh,$parent_id) = @_;
	
	my $sql = "select min(start) from cds where cds_parent_id = \'$parent_id\';";
	my $array_ref = $dbh->selectall_arrayref($sql);
	return $array_ref->[0]->[0];
};

sub get_cds_end_by_parent_id{
	my($dbh,$parent_id) = @_;
	
	my $sql = "select max(end) from cds where cds_parent_id = \'$parent_id\';";
	my $array_ref = $dbh->selectall_arrayref($sql);
	return $array_ref->[0]->[0];
};

sub get_cds_pos_by_parent_id{
	my($dbh,$parent_id) = @_;
	
	my $start = get_cds_start_by_parent_id($dbh,$parent_id);
	my $end   = get_cds_end_by_parent_id($dbh,$parent_id);
	
	return ($start,$end);
};
1;