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

Initialize

=head1 SYNOPSIS

	use Initialize;
	
	Initialize::load_gene_set($dbh,$config,$validation_file,$source,$gff_file,$fasta_file,$dsn,$user,$pass);

=head1 DESCRIPTION

This module loads the GFF into a SeqFeature Store database either in memory or disk if connection details is given.
The genes that passes the validation is inserted into the gene model table.

=head1 Author

	Mikkel B Christensen

=head1 METHODS

=cut

package Initialize;
use strict;
use warnings;
use v5.26;

use DBI;
use Data::Dumper;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);
use Log::Log4perl;
use Digest::MD5 qw(md5_hex);
use Validate;
use GeneModel qw(%BIOTYPE);

use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

my $verbose = 0;
my $db;
my $loader;

=head2 load_gene_set

 Title: load_gene_set
 Usage: Initialize::load_gene_set($dbh,$config,$validation_file,$source,$gff_file,$fasta_file,$dsn,$user,$pass)
 Function: load the gff into a SeqFeature Store database, then loads validated genes into the gene_model table.
 Returns: Hashref of stats for the number of genes, obsolete, etc.
 Args: Database handle object, config object,source of GFF, GFF file,FASTA file,database connection details (dns),user,pass
=cut

sub load_gene_set {
  my (
    $dbh,        $config, $validation_file, $source, $gff_file,
    $fasta_file, $dsn,    $user,            $pass,   $prot_fasta
  ) = @_;
  my %stats = (
    obsolete           => [],
    not_finished       => [],
    not_validated      => [],
    no_gene_model      => [],
    total_loaded       => 0,
    ok_gene_model      => 0,
    preloaded_features => 0,
    genes              => 0,
  );
  my $infoLog = Log::Log4perl->get_logger("infologger");

  my $proteins = _load_proteins($prot_fasta);

  $stats{preloaded_features} = _gff_load($gff_file, $fasta_file, $dsn, $user, $pass);

  _check_biotypes($db);

  my @genes      = $db->get_features_by_type('gene');
  my @prot_genes = $db->get_features_by_type('protein_coding_gene');
  my @pseudogenes = $db->get_features_by_type('pseudogene');
  my $fingerprints = {};

  $stats{genes} = scalar(@genes);
  @genes = (@genes, @prot_genes, @pseudogenes);
  open my $validation_fh, '>>', $validation_file;
  foreach my $gene (@genes) {
    my %attb    = $gene->attributes;
    my $gene_id = $attb{load_id}->[0];

    # Exclude some gene models
    if ($source eq 'cap') {
      if (_gene_is_obsolete(\%attb)) {
        $infoLog->info("Gene $gene_id is obsolete");
        push @{$stats{obsolete}}, $gene_id;
        next;
      }
      if (not _gene_is_finished(\%attb)) {
        $infoLog->info("Gene $gene_id is not finished");
        push @{$stats{not_finished}}, $gene_id;
        next;
      }
    }
    my $passed_validation = 1;
    $passed_validation = Validate::validate_gene($gene, $config, $validation_fh, $proteins, $fingerprints);

    if ($source eq 'cap' and $passed_validation < 0) {
      push @{$stats{not_validated}}, $gene_id;
      $infoLog->info("Gene $gene_id did not pass validation");
      next;
    }

    my ($gene_model, $CDS_present) = _build_gene_model($gene);

    if ($gene_model and %{$gene_model}) {
      $stats{ok_gene_model}++;
      _insert_gene_model($gene_model, $dbh, $source);
      _insert_exon($gene_model, $dbh, $source);
      if ($CDS_present) {
        _insert_CDS($gene_model, $dbh, $source);
      }
    } else {
      push @{$stats{no_gene_model}}, $gene_id;
    }

    $stats{total_loaded}++;
  }

  return \%stats;
}

sub _check_biotypes {
  my ($db) = @_;

  my %known_types = ();
  for my $name (keys %BIOTYPE) {
    %known_types = (%known_types, %{$BIOTYPE{$name}});
  }
  my @unsupported = ();
  for my $biotype ($db->types) {
    $biotype =~ s/:.+$//;
    if (not exists $known_types{$biotype}) {
      push @unsupported, $biotype;
    }
  }
  die("Unknown biotypes: " . join(", ", @unsupported)) if @unsupported;
}

sub _load_proteins {
  my ($fasta_path) = @_;
  return {} if not $fasta_path or not -s $fasta_path;

  say("Loading proteins...");
  my %prots;
  my $n       = 0;
  my $prot_in = Bio::SeqIO->new(-file => $fasta_path, -format => "Fasta");
  while (my $prot = $prot_in->next_seq()) {
    my $id = $prot->id;
    $prots{$id} = $prot->seq;
  }
  return \%prots;
}

sub _gene_is_obsolete {
  my ($attrib) = @_;
  my $obsolete = $attrib->{obsolete}->[0];

  return (defined($obsolete) and $obsolete eq 'true');
}

sub _gene_is_finished {
  my ($attrib) = @_;

  my $status = $attrib->{status}->[0];
  return (defined($status) and $status =~ /^Finished|Finished annotation$/);
}

sub _gff_load {
  my ($gff_file_name, $fasta_file_name, $dsn, $user, $pass) = @_;
  warn "load GFF $gff_file_name\n";
  $db     = ();
  $loader = ();
  if ($dsn) {
    $db = Bio::DB::SeqFeature::Store->new(
      -adaptor => 'DBI::mysql',
      -dsn     => $dsn,
      -user    => $user,
      -pass    => $pass,
      -create  => 1
    );
  } else {
    $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory',);
  }

  $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(
    -store   => $db,
    -verbose => $verbose,
  );

  open my $gff_fh, '<', $gff_file_name;
  my $loaded;
  if ($fasta_file_name) {
    open my $fasta_fh, '<', $fasta_file_name;
    $loaded = $loader->load($gff_fh, $fasta_fh);
  } else {
    my $loaded = $loader->load($gff_fh);
  }

  return $loaded;
}

sub _build_gene_model {
  my ($gene) = @_;
  my %attb = $gene->attributes;
  my %gene_model;
  my $CDS_present = 0;
  $gene_model{gene_id}               = $attb{load_id}->[0];
  $gene_model{validation_error_code} = $attb{validation_error_code}->[0];
  $gene_model{biotype} = $gene->type;
  $gene_model{biotype} =~ s/:.+$//;

  my @mRNAs = $gene->get_SeqFeatures('mRNA');
  my @transcripts = $gene->get_SeqFeatures();
  my %tr_biotypes = ();

  TRAN: foreach my $tr (@transcripts) {
    my %tr_model;
    my %rna_attb                   = $tr->attributes;
    my $RNA_validation_error_code  = $rna_attb{validation_error_code}->[0];

    $tr_model{transcript_id}         = $rna_attb{load_id}->[0];
    $tr_model{validation_error_code} = $RNA_validation_error_code;
    $tr_model{biotype} = $tr->type;
    $tr_model{biotype} =~ s/:.+$//;
    $tr_biotypes{$tr_model{biotype}} = 1;

    my @CDS = $tr->get_SeqFeatures('CDS');
    if (scalar(@CDS) > 0) {
      $CDS_present = 1;
      my $CDS_sequence = get_CDS($tr, '');
      return if not $CDS_sequence;

      my $prot_seq     = get_translation($CDS_sequence, $tr->strand);
      my $md5_checksum = md5_hex($prot_seq);

      my %cds_attb = $CDS[0]->attributes;
      $tr_model{cds_start}        = $CDS[0]->start;
      $tr_model{cds_end}          = $CDS[0]->end;
      $tr_model{CDS_Parent_id}    = $cds_attb{parent_id}->[0];
      $tr_model{CDS_md5_checksum} = $md5_checksum;
      $tr_model{cds_error_code}   = $RNA_validation_error_code;
    }

    my @exons = $tr->get_SeqFeatures('exon');

    if (@exons) {
      my $exon_number = 1;
      foreach my $exon (@exons) {
        my %exon_attb = $exon->attributes;
        my %exon_model;

        # Make unique exon id
        if ($exon_attb{load_id}->[0]) {
          $exon_model{exon_id} = $exon_attb{load_id}->[0];
        } else {
          $exon_model{exon_id} = $exon_attb{parent_id}->[0] . $exon_number;
          $exon_number++;
        }

        $exon_model{scaffold} = $exon->seq_id;
        $exon_model{strand}   = $exon->strand;
        $exon_model{start}    = $exon->start;
        $exon_model{end}      = $exon->end;

        push @{$tr_model{exon}}, \%exon_model;
      }
    } else {
      # No exons? Create one from the transcript
      my %exon_model = ();
      $exon_model{exon_id} = $tr_model{transcript_id} . "-exon";
      $exon_model{scaffold} = $tr->seq_id;
      $exon_model{strand}   = $tr->strand;
      $exon_model{start}    = $tr->start;
      $exon_model{end}      = $tr->end;

      push @{$tr_model{exon}}, \%exon_model;
    }
    push @{$gene_model{transcript}}, \%tr_model;
  }

  if ($gene_model{biotype} eq 'gene') {
    my $gene_biotype = join(",",(sort keys %tr_biotypes));
    $gene_model{biotype} = $gene_biotype;
  }

  return (\%gene_model, $CDS_present);
}

sub get_CDS {
  my ($mrna_obj, $validation_fh) = @_;
  my $infoLog = Log::Log4perl->get_logger("infologger");
  my $seq_id  = $mrna_obj->seq_id;
  my $strand  = $mrna_obj->strand;

  my @CDS  = $mrna_obj->get_SeqFeatures('CDS');
  return if not @CDS;
  my %attb = $CDS[0]->attributes;

  #only expects one CDS per mRNA, code will break if more.
  if (scalar @CDS > 1) {
    warn(
"more than one cds for mrna, only expects one cds per mrna, code will break if more. $attb{load_id}->[0]"
    );
    return;
  }
  my $CDS_start = $CDS[0]->start;
  my $CDS_end   = $CDS[0]->end;
  my $CDS_sequence;

  my @exons = $mrna_obj->get_SeqFeatures('exon');
  @exons = sort { $a->start <=> $b->start } @exons;

  #need to splice the mRNA
  if (scalar @exons == 1) {
    $CDS_sequence = get_sequence($seq_id, $CDS_start, $CDS_end);
  } elsif (scalar @exons > 1) {
    my $exon_start        = shift @exons;
    my $exon_end          = pop @exons;
    my $exon_start_endpos = $exon_start->end;
    my $exon_end_startpos = $exon_end->start;
    my $cds_in_one_exon   = 0;

    #do not include 5' UTR
    while ($exon_start_endpos < $CDS_start) {
      if (!scalar @exons) {
        $cds_in_one_exon = 1;
        last;
      }
      $exon_start        = shift @exons;
      $exon_start_endpos = $exon_start->end;
    }

    #do not include 3' UTR
    while ($CDS_end < $exon_end_startpos) {
      if (!scalar @exons) {
        $cds_in_one_exon = 1;
        last;
      }
      $exon_end          = pop @exons;
      $exon_end_startpos = $exon_end->start;
    }

    if ($cds_in_one_exon) {
      $CDS_sequence = get_sequence($seq_id, $CDS_start, $CDS_end);
    } else {
      my $CDS_start_seq = get_sequence($seq_id, $CDS_start, $exon_start_endpos);

      my $CDS_end_seq = get_sequence($seq_id, $exon_end_startpos, $CDS_end);

      chomp($CDS_start_seq, $CDS_end_seq);

      if (scalar(@exons)) {
        my $exon_sequence;
        foreach my $exon (@exons) {
          my $exon_start = $exon->start;
          my $exon_end   = $exon->end;
          $exon_sequence .= get_sequence($seq_id, $exon_start, $exon_end);

          #print $validation_fh "middle exon\t" . $exon_sequence . "\n";
        }
        chomp($CDS_start_seq, $exon_sequence, $CDS_end_seq);
        $CDS_sequence = $CDS_start_seq . $exon_sequence . $CDS_end_seq;
      } else {
        $CDS_sequence = $CDS_start_seq . $CDS_end_seq;
      }
    }
  }
  unless ($CDS_sequence) {
    $infoLog->info("NO CDS sequence was returned for mRNA, $attb{load_id}->[0]");
  }
  return $CDS_sequence;
}

sub get_translation {
  my ($CDS_sequence, $strand) = @_;

  return if not $CDS_sequence;

  my $seqobj = Bio::Seq->new(-seq => $CDS_sequence);
  if ($strand eq '-' or $strand eq '-1') {
    $seqobj = $seqobj->revcom();
  }
  my $prot_seq = $seqobj->translate->seq;

  return $prot_seq;
}

sub get_sequence {
  my ($scaffold, $start, $end) = @_;
  my $seq = $db->fetch_sequence($scaffold, $start, $end);
  return $seq;
}

sub _insert_exon {
  my ($hash, $dbh, $source) = @_;

  my $insert_sql = "insert exon(exon_id,
								  source,
								  scaffold,
								  strand,
								  start,
								  end
					 )
					 select ?,?,?,?,?,?;";
  my $sth = $dbh->prepare($insert_sql);
  foreach my $transcript (@{$hash->{transcript}}) {
    foreach my $exon (@{$transcript->{exon}}) {
      $sth->bind_param(1, $exon->{exon_id});
      $sth->bind_param(2, $source);
      $sth->bind_param(3, $exon->{scaffold});
      $sth->bind_param(4, $exon->{strand});
      $sth->bind_param(5, $exon->{start});
      $sth->bind_param(6, $exon->{end});

      $sth->execute();
    }
  }

}

sub _insert_CDS {
  my ($hash, $dbh, $source) = @_;

  my $insert_sql = "insert cds(cds_parent_id,
								start,
								end,
								md5_checksum,
								cds_error_code
					 )
					 select ?,?,?,?,?;";

  my $sth = $dbh->prepare($insert_sql);
  foreach my $transcript (@{$hash->{transcript}}) {

    $sth->bind_param(1, $transcript->{CDS_Parent_id});
    $sth->bind_param(2, $transcript->{cds_start});
    $sth->bind_param(3, $transcript->{cds_end});
    $sth->bind_param(4, $transcript->{CDS_md5_checksum});
    $sth->bind_param(5, $transcript->{cds_error_code});

    $sth->execute();
  }
}

sub _insert_gene_model {
  my ($hash, $dbh, $source) = @_;
  my $insert_sql = "insert gene_model(
										exon_id,
										transcript_id,
										gene_id,
										source,
										error_code,
                    biotype
									   )
							select ?,?,?,?,?,?;";
  my $sth     = $dbh->prepare($insert_sql);
  my $gene_id = $hash->{gene_id};
  foreach my $transcript (@{$hash->{transcript}}) {
    foreach my $exon (@{$transcript->{exon}}) {
      $sth->bind_param(1, $exon->{exon_id});
      $sth->bind_param(2, $transcript->{transcript_id});
      $sth->bind_param(3, $hash->{gene_id});
      $sth->bind_param(4, $source);
      $sth->bind_param(5, $hash->{validation_error_code});
      $sth->bind_param(6, $hash->{biotype});

      $sth->execute();
    }
  }
}

1;
