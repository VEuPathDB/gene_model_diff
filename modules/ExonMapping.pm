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

package ExonMapping;

use strict;
use warnings;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use Try::Tiny;
use Exon;
use lib('.');
use GeneModel;


sub work_out_exon_mapping{
	my ($dbh) = @_;
	my @exon_matches = qw(identical included partial_5 partial_3 spanning no_match_cap no_match_vb);#no_math_cap no_macth_vb
	my %mapped_exons;
	my $errorLog = Log::Log4perl->get_logger("ExonMapping::work_out_exon_mapping");
	foreach my $exon_match (@exon_matches){
		try{
			if($exon_match eq 'no_match_vb'){
				find_unmatch_vb_exon($dbh);			
			}else{
				my $over_lap_ref = Exon::get_overlapping_exons($dbh,$exon_match);
				insert_exon_mappings($dbh,\%mapped_exons,$over_lap_ref,$exon_match);
			}
		}catch{
			$errorLog->error($_);	
		}
			

	}	
}

sub insert_exon_mappings{
	my ($dbh,$mapped_exons,$array_ref,$macth) = @_;
	my $insert_sth = $dbh->prepare(get_sql('insert_exon_mapping'));
	my $size = keys %{$mapped_exons};
	
	
	foreach my $row (@{$array_ref}){	
		if(!$row->[1]){
			$row->[1] = 'none';
		}
		my $key = "$row->[0]_$row->[1]";
		if((exists $mapped_exons->{$key})){			
			next;
		}else{
			$mapped_exons->{$key} = 1;	
		}
			
		if($row->[0] and $row->[0] ne 'vb'){
			$insert_sth->bind_param(1,$row->[0]);	
		}
		
		if($row->[1] and $row->[1] ne 'cap'){			
			$insert_sth->bind_param(2,$row->[1]);			
		}
		
		$insert_sth->bind_param(3,$macth);		
		$insert_sth->execute();		
	}
	
}

sub find_unmatch_vb_exon {
	my($dbh) = @_;
	my $mapped_vb_exons_sth = $dbh->prepare(get_sql('get_mapped_vb_exons'));
	my $mapped_vb_transcript_sth = $dbh->prepare(get_sql('get_mapped_vb_transcript'));
	
	my $insert_non_mapped_exons_sth = $dbh->prepare(get_sql('insert_non_mapped_exons'));
	
	$mapped_vb_exons_sth->execute();
	
	my $mapped_vb_exons_array_ref = $mapped_vb_exons_sth->fetchall_arrayref;
	
	foreach my $row (@{$mapped_vb_exons_array_ref}){		
		
		my $mapped_vb_transcript_id = GeneModel::get_gene_model_by_id($dbh,'exon',$row->[0],'vb','transcript');
				
		my $vb_exons_of_mapped_vb_transcript_array_ref =  GeneModel::get_gene_model_by_id($dbh,'transcript',$mapped_vb_transcript_id,'vb');  
		foreach my $row (@{$vb_exons_of_mapped_vb_transcript_array_ref}){				
				$insert_non_mapped_exons_sth->execute($row->[0],$row->[0]);
		}
				
	}	
}

sub get_exon_mappings_by_id {
	my($dbh,$exon_type,$id) = @_;
	
	my %exon_types =( "cap" => "cap_exon_id",
					  "vb"  => "vb_exon_id"
					  );
	
	my $sql = "select cap_exon_id,vb_exon_id,map_type from exon_mappings, gene_model where $exon_types{$exon_type} = exon_id and transcript_id =\'$id\';";
	my $array_ref = _submit_sql($dbh,$sql);	
	
	return $array_ref;
}

sub get_exon_mappings {
	my($dbh,$exon_type,$id) = @_;
		
	my %exon_types =( "cap" => "cap_exon_id",
					  "vb"  => "vb_exon_id"
					  );
	
	my $return_type;
	if($exon_type eq 'cap'){
		$return_type = 'vb';
	}else{$return_type = 'cap'}
	
	my $sql;
	if($exon_type and $id){
		$sql = "select map_type, $exon_types{$return_type} from exon_mappings where $exon_types{$exon_type} = \'$id\'";
	}else{
		$sql = "select cap_exon_id,vb_exon_id,map_type from exon_mappings;";
	}
	#warn $sql;
	my $array_ref = _submit_sql($dbh,$sql);
	return $array_ref;
	
}

sub _submit_sql {
	my($dbh,$sql) = @_;
	my $array_ref = $dbh->selectall_arrayref($sql);
	return $array_ref;	
}

sub get_sql{
	my ($sql_name) = @_;
	
	if($sql_name eq 'insert_exon_mapping'){
			my $insert_mapping_sql = "insert exon_mappings(
									   cap_exon_id,
									   vb_exon_id,
									   map_type)
									   select ?,?,?;";
									   
					 return $insert_mapping_sql;
	
	}elsif($sql_name eq 'insert_non_mapped_exons'){
			my $non_mapped_exons = "insert exon_mappings(cap_exon_id,vb_exon_id,map_type) select 'none',exon_id,'no_match_vb' from gene_model where exon_id = ? and not exists (select * from exon_mappings where vb_exon_id = ?);";
			return $non_mapped_exons;
	}elsif($sql_name eq 'get_mapped_vb_exons'){
		my $mapped_vb_exons = "select vb_exon_id from exon_mappings where vb_exon_id is not NULL;";
			return $mapped_vb_exons;		
	}
}

1;