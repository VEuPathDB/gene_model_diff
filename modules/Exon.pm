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

package Exon;

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;


=head2 get_overlapping_exons

 Title:    get_overlapping_exons	
 Usage:    Exon::get_overlapping_exons().	
 Function: Execute SQL on the database. 	
 Returns:  Array refference of pairs of overlapping exon ids. 	
 Args:     Database handle object,type of overlap, where type can be one of the following strings (identical|partial_3|partial_5|spanning|included|no_match_cap). 
=cut 

sub get_overlapping_exons{
	my ($dbh,$overlap_type) = @_;
	my $sql = get_sql($overlap_type);
	#log error if no SQL
	unless($sql){
		confess("There is no SQL of type $overlap_type, therefore no Exons will be returned");
	}
	my $sth = $dbh->prepare($sql);
	
	$sth->execute();
	my $array_ref = $sth->fetchall_arrayref;
	
	return $array_ref;	
}



sub get_sql{
	my ($sql_name) = @_;
	
	if($sql_name eq 'identical'){
		my $identical_sql = "select cap.exon_id, vb.exon_id from exon cap, exon vb
						 where cap.source   = 'cap'
						 and   vb.source    = 'vb'
						 and   cap.scaffold = vb.scaffold
						 and   cap.strand   = vb.strand
						 and   cap.start    = vb.start
						 and   cap.end      = vb.end;";
						 
						 return $identical_sql;						 
	}elsif($sql_name eq 'partial_3'){				 
		my $partial_3_sql = "select cap.exon_id, vb.exon_id from exon cap, exon vb
					 	where cap.source   = 'cap'
					 	and   vb.source    = 'vb'
					 	and   cap.scaffold = vb.scaffold
					 	and   cap.strand   = vb.strand
					 	and   cap.end   >= vb.start
					 	and   cap.end   <= vb.end;";

                         return $partial_3_sql						 						 
	}elsif($sql_name eq 'partial_5'){
		my $partial_5_sql = "select cap.exon_id, vb.exon_id from exon cap, exon vb
						 where cap.source   = 'cap'
						 and   vb.source    = 'vb'
						 and   cap.scaffold = vb.scaffold
						 and   cap.strand   = vb.strand
						 and   cap.start >= vb.start
						 and   cap.start <= vb.end;";
						 
						 return $partial_5_sql;
					 
	}elsif($sql_name eq 'spanning'){				 
		my $spanning_sql = "select cap.exon_id, vb.exon_id from exon cap, exon vb
					 where cap.source   = 'cap'
					 and   vb.source    = 'vb'
					 and   cap.scaffold = vb.scaffold
					 and   cap.strand   = vb.strand
					 and   cap.start   <= vb.start
					 and   cap.end     >= vb.end;";
					 
					 return $spanning_sql;
					 
	}elsif($sql_name eq 'included'){				 
		my $included_sql = "select cap.exon_id, vb.exon_id from exon cap, exon vb
					 where cap.source   = 'cap'
					 and   vb.source    = 'vb'
					 and   cap.scaffold = vb.scaffold
					 and   cap.strand   = vb.strand
					 and   cap.start   >= vb.start
					 and   cap.end     <= vb.end;";
					 
					 return $included_sql;
					 
	}elsif($sql_name eq 'no_match_cap'){
		my $no_match_cap_sql = "select exon_id,source from exon where source = 'cap' and  not exists (select * from exon_mappings where exon_id = cap_exon_id);";
		return $no_match_cap_sql;
	}
}
1;