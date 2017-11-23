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

package GeneModel;
use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use DBI;

my %types = ("gene"       => "gene_id",
			 "transcript" => "transcript_id",
			 "exon"       => "exon_id"
			);

my %return_types = ("gene"        => "2",
			 		"transcript" => "1",
			 		"exon"       => "0",
			 		"gene_model" => "0..2"
				   );

sub get_gene_model_by_id {
	my($dbh,$type,$id,$source,$return_type,$force_array) = @_;
	confess("missing parameter to get_gene_model_by_id") unless($dbh and $type and $id);
	my $sql;
	if($source){
		$sql = "select exon_id,transcript_id,gene_id,source
			   	   from gene_model where $types{$type} = \'$id\' and source = \'$source\';";
	}else{	
		$sql = "select exon_id,transcript_id,gene_id,source
			   	   from gene_model where $types{$type} = \'$id\';";
   }
   
   my $array_ref = _submit_sql($dbh,$sql);
   
  
   
   if(scalar @{$array_ref} == 1 and $return_type and $force_array){
   	   my @IDs;
   	   foreach my $row (@{$array_ref}){
   	   	   my $id = @{$row}[$return_types{$return_type}];
   	   	   push @IDs,$id;
   	   }
   	   return \@IDs;  
   }elsif(scalar @{$array_ref} == 1 and $return_type){ 
   	   return $array_ref->[0][$return_types{$return_type}];
   	   
   }elsif(scalar @{$array_ref} > 1 and $return_type){
   	   my @IDs;
   	   foreach my $row (@{$array_ref}){
   	   	   my $id = @{$row}[$return_types{$return_type}];
   	   	   push @IDs,$id;
   	   }
   	   return \@IDs;  	   
   }else{
   	   return $array_ref;
   }	
}

sub get_distinct_id_by_id {
	my($dbh,$type,$id,$return_id,$source,$force_array) = @_;
	confess("missing parameter to get_distinct_id_by_id") unless($dbh and $type and $id and $return_id and $source);
	my $sql = "select distinct $types{$return_id} from gene_model where $types{$type} = \'$id\' and source = \'$source\';";
	my $array_ref = _submit_sql($dbh,$sql);	
	
	if($force_array){
		return $array_ref;
	}else{return $array_ref->[0]->[0];}
	
}

sub get_distinct_id_by_source {
	my($dbh,$type,$source) = @_;
	
	my $sql = "select distinct $types{$type} from gene_model where source = \'$source\';";
	my $array_ref = _submit_sql($dbh,$sql);	
	return $array_ref;
}

sub _submit_sql {
	my ($dbh,$sql) = @_;
	my $array_ref = $dbh->selectall_arrayref($sql);
	return $array_ref;	
}

1;