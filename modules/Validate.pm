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

Validate

=head1 SYNOPSIS

	use Validate;
	
	Validate::validate_gene($gene,$config,$validation_fh);

=head1 DESCRIPTION

This module validates the gene model. The following checks are applied: All sub-feature
must be on the same strand and within its parents boundaries. Coding sequence must have Start and Stop codon
or have a 'no-ATG'/'no-stop' tag to pass.

=head1 Author

	Mikkel B Christensen

=head1 METHODS

=cut

package Validate;

use strict;
use warnings;
use v5.26;
use autodie qw(:all);
use Carp qw(cluck carp croak confess);

use Log::Log4perl;
use Initialize;

use Bio::Seq;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

=head2 validate_gene

 Title: validate_gene
 Usage: Validate::validate_gene($gene,$config,$validation_fh)
 Function: validates the gene model
 Returns: validation/error code, which are given in the config file. 
 Args: SeqFeature gene object, config object, validation file handle
=cut

sub validate_gene {
  my ($gene, $config, $validation_fh, $proteins, $cds_fingerprints) = @_;
  my %gene_attb = $gene->attributes;
  my $errorLog  = Log::Log4perl->get_logger("errorlogger");

  my %validation_error_code = (
    'GFF_mRNA'          => $config->val('Validate', 'GFF_mRNA'),
    'GFF_exon'          => $config->val('Validate', 'GFF_exon'),
    'CDS_start'         => $config->val('Validate', 'CDS_start'),
    'CDS_stop'          => $config->val('Validate', 'CDS_stop'),
    'CDS_internal_stop' => $config->val('Validate', 'CDS_internal_stop'),
    'runtime_error'     => $config->val('Validate', 'runtime_error')
  );

  my @approved_users            = $config->val('Validate', 'approved_user');
  my %validation_approved_email = map { $_ => 'approved' } @approved_users;

  my %gene_hash;
  $gene_hash{scaffold} = $gene->seq_id;
  $gene_hash{strand}   = $gene->strand;
  $gene_hash{start}    = $gene->start;
  $gene_hash{end}      = $gene->end;
  my $strand = $gene->strand;
  my $gene_id = $gene_attb{load_id}->[0];

  my @RNAs                   = $gene->get_SeqFeatures('mRNA');
  my $gene_validation_status = 0;
  my $validation_string      = '';
  my $exon_fingerprints = {};
  my @gene_cds_fingerprints;

  foreach my $mRNA (@RNAs) {
    my %mRNA_hash;
    my %mRNA_attb              = $mRNA->attributes;
    my $mRNA_id                = $mRNA_attb{load_id}->[0];
    my $prot_fasta_seq         = $proteins->{$mRNA_id} ? $proteins->{$mRNA_id} . "*" : "";
    my $mRNA_validation_status = 0;

    my $NO_ATG        = 0;
    my $NO_STOP       = 0;
    my $INTERNAL_STOP = 0;

    # Check that the mRNA has no more than 1 CDS
    my $cds_ids_count = _check_cds_count($mRNA);
    if ($cds_ids_count > 1) {
        $validation_string .= "mRNA:$mRNA_id|Too many CDS ids: $cds_ids_count;";
        $mRNA_validation_status += $validation_error_code{GFF_mRNA};
    }

    # Check CDS fingerprints (within all genes)
    push @gene_cds_fingerprints, _get_mRNA_fingerprint($mRNA, ['CDS']);

    # Check exons and CDS fingerprints (only within gene)
    my $exon_fingerprint = _get_mRNA_fingerprint($mRNA, ['CDS', 'exon']);
    if ($exon_fingerprint and exists $exon_fingerprints->{$exon_fingerprint}) {
        $validation_string .= "mRNA:$mRNA_id|Duplicate_mRNA;";
        $mRNA_validation_status += $validation_error_code{GFF_mRNA};
    } else {
        $exon_fingerprints->{$exon_fingerprint} = 1;
    }

    my $has_no_atg  = exists $mRNA_attb{'no-ATG'};
    my $has_no_stop = exists $mRNA_attb{'no-STOP'};
    my $has_author =
      (exists $mRNA_attb{owner} and $validation_approved_email{$mRNA_attb{owner}->[0]});
    my $is_left_partial =
      (exists $mRNA_attb{'is_fmin_partial'} and ($mRNA_attb{'is_fmin_partial'}->[0] eq 'true'));
    my $is_right_partial =
      (exists $mRNA_attb{'is_fmax_partial'} and ($mRNA_attb{'is_fmax_partial'}->[0] eq 'true'));
    my $is_start_partial = (
           ($strand == 1 and $is_left_partial)
        or ($strand == '-1' and $is_right_partial)
    );
    my $is_end_partial = (
           ($strand == 1 and $is_right_partial)
        or ($strand == '-1' and $is_left_partial)
    );
    my $stop_seleno = (exists $gene_attb{'Note'}
        and ($gene_attb{'Note'}->[0] eq 'stop_codon_redefined_as_selenocysteine'));

    # Check both missing start/stop codon
    # and allow if partial, or owner approved
    if ($has_no_atg) {
      if ($has_author) {
        $NO_ATG = 2;
      } else {
        $NO_ATG = 1;
      }
    }
    if ($is_start_partial) {
      $NO_ATG = 2;
    }

    if ($has_no_stop) {
      if ($has_author or $is_end_partial) {
        $NO_STOP = 2;
      } else {
        $NO_STOP = 1;
      }
    }
    if ($is_end_partial) {
      $NO_STOP = 2;
    }
    if ($stop_seleno) {
      $INTERNAL_STOP = 2;
    }

    eval {
      my $mRNA_ID = $mRNA_attb{load_id}->[0];
      $mRNA_hash{scaffold} = $mRNA->seq_id;
      $mRNA_hash{strand}   = $mRNA->strand;
      $mRNA_hash{start}    = $mRNA->start;
      $mRNA_hash{end}      = $mRNA->end;

      my $validation_text;

      # Check mRNA errors
      if (!_check_subfeature(\$validation_text, $mRNA_ID, \%gene_hash, \%mRNA_hash, 'Gene:mRNA')) {
        $validation_string .= "mRNA:$mRNA_ID|$validation_text;";
        $mRNA_validation_status += $validation_error_code{GFF_mRNA};
      }

      # Check exons errors
      my @exons      = $mRNA->get_SeqFeatures('exon');
      my $exon_error = 0;
      foreach my $exon (@exons) {
        my %exon_hash;
        my %exon_attb = $exon->attributes;
        my $exon_ID   = $exon_attb{load_id}->[0];
        $exon_hash{scaffold} = $exon->seq_id;
        $exon_hash{strand}   = $exon->strand;
        $exon_hash{start}    = $exon->start;
        $exon_hash{end}      = $exon->end;

        if (!_check_subfeature(\$validation_text, $exon_ID, \%mRNA_hash, \%exon_hash, 'mRNA:exon'))
        {
          $validation_string .= "exon:$exon_ID|$validation_text;";
          $exon_error = 1;
        }

      }

      if ($exon_error) {
        $mRNA_validation_status += $validation_error_code{GFF_exon};
      }

      # Check CDS sequence errors
      my $CDS_sequence = Initialize::get_CDS($mRNA, $validation_fh);
      my $prot_seq = "";
      unless ($CDS_sequence) {
        $errorLog->error("No CDS for $mRNA_ID");
      } else {
        $prot_seq = Initialize::get_translation($CDS_sequence, $mRNA_hash{strand});
        unless ($prot_seq) {
          $errorLog->error("No translation for $mRNA_ID");
        }
      }

      my $sequence_check = "";

      if (not $prot_seq) {
        $sequence_check="no_sequence";
        $validation_string .= "No_sequence;";
        $mRNA_validation_status += $validation_error_code{GFF_mRNA};
      } elsif (not $prot_fasta_seq or $prot_seq eq $prot_fasta_seq) {
        $sequence_check = _check_seq($prot_seq, \$validation_text);
      } else {
        #say("Protein sequence differs:\n>$mRNA_id\n$prot_seq\n----\n$prot_fasta_seq");
        $sequence_check = _check_seq($prot_fasta_seq, \$validation_text);
      }

      if ($sequence_check ne 'passed' and $sequence_check ne 'no_sequence') {
        $validation_string .= "mRNA:$mRNA_ID";
        $validation_string .= "$validation_text;";
        $errorLog->error($validation_string);
        my ($start_flag, $stop_flag, $internal_stop_count) =
          split /:/, $sequence_check;
        if (!$start_flag) {
          $mRNA_validation_status += $validation_error_code{CDS_start};
        }

        if (!$stop_flag) {
          $mRNA_validation_status += $validation_error_code{CDS_stop};
        }

        if (!$internal_stop_count == 0) {
          $mRNA_validation_status += $validation_error_code{CDS_internal_stop};
        }

      }

      $mRNA->add_tag_value('validation_error_code' => $mRNA_validation_status);
      $mRNA->update();
    };

    #warn "$mRNA_validation_status : $validation_error_code{CDS_stop} : $NO_STOP";
    if ($@) {
      $errorLog->error($@);
      say("Check failed: $@");
      $gene_validation_status = -1;

      # We have some mRNA errors to check
    } elsif ($mRNA_validation_status
      and $gene_validation_status > -1
      and $gene_validation_status < $mRNA_validation_status)
    {

      # Approve missing CDS start codon
      if (  $mRNA_validation_status == $validation_error_code{CDS_start}
        and $NO_ATG == 2)
      {
        $gene_validation_status = 0;

        # Approve missing stop codon
      } elsif ($mRNA_validation_status == $validation_error_code{CDS_stop}
        and $NO_STOP == 2)
      {
        $gene_validation_status = 0;

        # Approve missing stop codon
      } elsif ($mRNA_validation_status == $validation_error_code{CDS_internal_stop}
        and $INTERNAL_STOP == 2)
      {
        $gene_validation_status = 0;

        # Approve both missing stop and start codon
      } elsif (
        (
          $mRNA_validation_status ==
          ($validation_error_code{CDS_stop} + $validation_error_code{CDS_start})
        )
        and $NO_ATG == 2
        and $NO_STOP == 2
        )
      {
        $gene_validation_status = 0;

        # Approve both internal stop and start codon
      } elsif (
        (
          $mRNA_validation_status ==
          ($validation_error_code{CDS_internal_stop} + $validation_error_code{CDS_start})
        )
        and $NO_ATG == 2
        and $INTERNAL_STOP == 2
        )
      {
        $gene_validation_status = 0;

        # Add the mRNA error to the gene
        # Missing CDS start codon
      } elsif ($mRNA_validation_status == $validation_error_code{CDS_start}
        and $NO_ATG == 1)
      {
        $gene_validation_status = $mRNA_validation_status;

        # Missing stop codon
      } elsif ($mRNA_validation_status == $validation_error_code{CDS_stop}
        and $NO_STOP == 1)
      {
        $gene_validation_status = $mRNA_validation_status;

        # Missing both start and stop codon
      } elsif (
        (
          $mRNA_validation_status ==
          ($validation_error_code{CDS_stop} + $validation_error_code{CDS_start})
        )
        and $NO_ATG == 1
        and $NO_STOP == 1
        )
      {
        $gene_validation_status = $mRNA_validation_status;

        # Other?
      } else {
        $gene_validation_status = -$mRNA_validation_status;
      }
    }
  }

  # Check CDS fingerprints
  my %gene_fingerprints = map { $_ => 1 } @gene_cds_fingerprints;
  for my $cds_fingerprint (sort keys %gene_fingerprints) {
    if ($cds_fingerprint and exists $cds_fingerprints->{$cds_fingerprint}) {
        $validation_string .= "Duplicate_gene_CDS=$cds_fingerprint;";
        $gene_validation_status = -$validation_error_code{GFF_mRNA};
    } else {
        $cds_fingerprints->{$cds_fingerprint} = 1;
    }
  }

  # Log if any error found
  if ($gene_validation_status) {
    my $gene_id = $gene_attb{load_id}->[0];
    my $owner   = $gene_attb{owner}->[0] || 'None';

    print $validation_fh "gene:$gene_id\t$owner\t$validation_string\n";
  }
  $gene->add_tag_value('validation_error_code' => $gene_validation_status);
  return $gene_validation_status;
}

sub _check_cds_count {
  # Return the number of distinct CDS ids for this mRNA
  my ($mRNA) = @_;

  my %cds_ids = ();
  for my $cds ($mRNA->get_SeqFeatures('CDS')) {
    my %cds_attb = $cds->attributes;
    my $cds_id = $cds_attb{load_id}->[0];
    $cds_ids{$cds_id} = 1;
  }

  return scalar %cds_ids;
}

sub _get_mRNA_fingerprint {
  my ($mRNA, $names) = @_;

  my @fingerprints;
  for my $name (@$names) {
    my @feats;
    for my $feat (sort $mRNA->get_SeqFeatures($name)) {
      push @feats, $feat->seq_id. ":" . $feat->strand . ":" . $feat->start . "-" . $feat->end;
    }
    my $feat_fingerprint = @feats ? "$name=" . join(",", sort @feats) : "";
    push @fingerprints, $feat_fingerprint if $feat_fingerprint;
  }

  return join(";", @fingerprints);
}

sub _check_subfeature {
  my ($validation_text, $ID, $feature, $sub_feature, $key) = @_;
  my $passed = 1;
  my ($feature_name, $sub_feature_name) = split /:/, $key;
  $$validation_text = '';
  if ($feature->{scaffold} ne $sub_feature->{scaffold}) {
    $$validation_text .=
"error:Scaffold not identical|proof:[$feature_name scaffold $feature->{scaffold}, $sub_feature_name scaffold $sub_feature->{scaffold}];";
    $passed = 0;
  }

  if ($feature->{strand} ne $sub_feature->{strand}) {
    $$validation_text .=
"error:Strand not identical|proof:[$feature_name strand $feature->{strand}, $sub_feature_name strand $sub_feature->{strand}];";
    $passed = 0;
  }

  if ($feature->{start} > $sub_feature->{start}) {
    $$validation_text .=
"error:sub feature start is not within borders|proof:[$feature_name start $feature->{start}, $sub_feature_name start $sub_feature->{start}];";
    $passed = 0;
  }

  if ($feature->{end} < $sub_feature->{end}) {
    $$validation_text .=
"error:sub feature end is not within borders|proof:[$feature_name end $feature->{end}, $sub_feature_name end $sub_feature->{end}];";
    $passed = 0;
  }
  return $passed;
}

sub _check_seq {
  my ($prot_seq, $validation_text) = @_;
  my $start_aa            = substr($prot_seq, 0,  1);
  my $end_aa              = substr($prot_seq, -1, 1);
  my $start_flag          = 0;
  my $stop_flag           = 0;
  my $internal_stop_count = 0;
  if ($start_aa eq 'M') { $start_flag++; }
  if ($end_aa eq '*')   { $stop_flag++; }
  $internal_stop_count -= $stop_flag;
  while ($prot_seq =~ m/\*/g) { $internal_stop_count++; }
  $$validation_text = '';

  if (!$start_flag) {
    $$validation_text .= "|error:No start codon";
  }
  if (!$stop_flag) {
    $$validation_text .= "|error:No stop codon";
  }
  if ($internal_stop_count > 0) {
    $$validation_text .= "|error:Internal stop codon";
  }

  if ($$validation_text) {
    $$validation_text .= "|proof:$prot_seq;";
  }

  if ($start_flag and $stop_flag and ($internal_stop_count == 0)) {
    return 'passed';
  } else {
    return "$start_flag:$stop_flag:$internal_stop_count";
  }
}
1;
