=pod 

=head1 NAME

Bio::EnsEMBL::Hive::PHI_runnables::EntryValidator

=head1 SYNOPSIS

    standaloneJob.pl Bio::EnsEMBL::Hive::PHI_runnables::EntryValidator  -debug 1

=head1 DESCRIPTION

    This runnable takes a validated entry from the phi_base csv file and queries the Taxonomy DB for the species node of the taxa identifier. 
    It returns also the associated children of the parent node (all substrains)

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2017] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    Please subscribe to the Hive mailing list:  http://listserver.ebi.ac.uk/mailman/listinfo/ehive-users  to discuss Hive-related questions or to be notified of our updates

=cut


package Bio::EnsEMBL::Hive::PHI_runnables::EntryValidator;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Scalar::Util qw(reftype);

sub param_defaults {
    return {'gzip_flags' => ""};
}

sub pre_cleanup {
    # This is where the runnable would take care of getting things
    # (typically files or database connections) in order before starting to run.
    # This is often useful if the jobs sometimes need to be retried before
    # completing successfully
}

sub fetch_input {
    my $self = shift;
    # read_input is typically used to validate input, test db connections are
    # working, etc. Here, we'll check that we were passed a parameter called
    # 'output_ids'

    
    my $annotn_tax_id = $self->param_required("_6");
    if (! $annotn_tax_id) {
         die 'annotn_tax_id is not defined'; # Will cause job to fail and leave a message in log_message
    }

    
    # use param_required for this, since we want the job to fail if there's no output_ids.
    # my $output_ids = $self->param_required('input_id');


    # if (! $output_ids) {
    #     die 'output_ids is not defined'; # Will cause job to fail and leave a message in log_message
    # }

}


sub run {
    my $self = shift;
    
    
    my $annotn_tax_id = $self->param_required("_6");

    ####### THIS HERE PROBABLY NEEDS TO BE DEREFERENCED TO SOMETHING!!!!!!!!!
    print "####### THIS HERE PROBABLY NEEDS TO BE DEREFERENCED TO SOMETHING!!!!!!!!! line 92 SubtaxonsFinder.pm\n";
    my $branch_dbas = $self->param_required("_15");

    my $search = Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();
    $self->param('search', $search);

    for my $brch_dba (@{$branch_dbas}) {

      my $mc = $brch_dba->get_MetaContainer();
      my $branch_species = $mc->single_value_by_key('species.production_name');
      my $division =  $mc->single_value_by_key('species.division');

      print "$branch_species \n";

    #   if ( $division eq $opts->{division} ) {
    #     # print "\t === Found genome for species " . $branch_species . " in " . $brch_dba->dbc()->dbname() . "/" . $brch_dba->species_id ."\n";
    #     # if ($is_first_brch_dba eq 'true') {
    #     #   print ('-FIRST DBA-');
    #     # } else {
    #     #   print('-SECONDARY DBA-');
    #     # }

    #     ( my $translation , my $id_type ) = $self->find_translation( $brch_dba, $annotn_uniprot_acc, $annotn_locus, $annotn_gene_name );
 
    #     if ($translation) {
    #       # print "\t\t!! Found gene on ensembl : ";
    #       my $gene_adaptor        = $brch_dba->get_adaptor("Gene");
    #       my $gene = $gene_adaptor->fetch_by_translation_stable_id($translation->stable_id);
          
    #       # print "using $id_type\n\t\t" . $gene->stable_id . "\t" . $translation->stable_id . "\n";
    #       # print $output $result_line . ',0,' . $id_type . ',' . $branch_species . "," . $translation->stable_id . ',' . 100 . "\n";

    #       next LINE;

    #     } else {
    #       # print "\t\tDidn't find direct match on Ensembl \n";
    #       my $up;
    #       eval {
    #         $up = $self->get_uniprot_seq($annotn_uniprot_acc);
    #         $uniprots->{$annotn_uniprot_acc} = $up;
    #         $id_type = 'uniprot';
    #       };
    #       if($@ || $uniprots->{$annotn_uniprot_acc}->{seq} eq '') {
    #         print "Could not find entry : " . $annotn_uniprot_acc;
    #         $up = $self->get_uniparc_seq($annotn_uniprot_acc);
    #         $uniprots->{$annotn_uniprot_acc} = $up;
    #         $id_type = 'uniparc';
    #       }

    #       my $count_seqs = $up->{seq}   =~ tr/>//; # insures that uniprot/uniparc queries return no more than one protein sequence. Skips the entry otherwise
    #       if ($count_seqs > 0) {
    #         print $error_output " *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- ERROR : skipping $annotn_phi_entry  *-*-*-*-*-*-*-*-*-*-*-*-*-*) \n";
    #         print $error_output " $id_type  query returns $count_seqs extra sequence(s): " . $up->{seq} . " \n";
    #         next LINE;
    #       }

    #       if ($uniprots->{$annotn_uniprot_acc}->{seq} ne '') { 

    #         print "\t-Got uniprot/uniparc accession.- \n";
    #         if ($is_first_brch_dba eq 'true' ) {
    #           print Dumper($uniprots->{$annotn_uniprot_acc});
    #         }

    #         # Create BLAST DB for this specie
    #         my $blast_db_dir = $opts->{blast_db_dir};
    #         my $blast_database = $blast_db_dir . "/" .$branch_species . ".pep.all";
    #         print "BLAST search: creating peptide file if it doesnt exist already in $blast_database:\n";
          
    #         if ( ! -f $blast_database ) { # Returns undef if the file doesn't exist, or if the file exists but it's not plain file (e.g. if it's a directory).
    #           open (my $db_file, ">", $blast_database) or croak "Could not open ".$blast_database;;
    #           print "getting all genes for $branch_species \n";

    #           my $translation_adaptor=$brch_dba->get_TranslationAdaptor();
    #           my @translations = @{$translation_adaptor->fetch_all()}; 
    
    #           foreach my $translation (@translations) {
    #             my $stable_id = $translation->stable_id();
    #             my $sequence = $translation->seq();
    #             print $db_file ">" . $stable_id . "\n";
    #             print $db_file $sequence . "\n";
    #           }
    #           close $db_file;
    #         } 

    #         # Create a BlastPlus factory
    #         my $fac;
           
    #         eval {
    #           $fac = Bio::Tools::Run::StandAloneBlastPlus->new( 
    #             -db_name => $branch_species,
    #             -db_dir => $blast_db_dir ,
    #             -db_data => $blast_database,
    #             -create => 1
    #           );
    #           $fac->make_db();
    #         };
    #         $@ =~ /EXCEPTION|WARNING/ and my $e = $@;
    #         if ( defined $e ) {
    #           print "\t******************* WARNING: **************************\n";
    #           print "BlastPlus factory could not be created for $annotn_uniprot_acc: This entry will be skipped\n";
    #           next LINE;
    #         }

    #         my $query   = Bio::Seq->new( -display_id => $annotn_uniprot_acc, -seq =>  $uniprots->{$annotn_uniprot_acc}->{seq} );
    #         my $results = $fac->blastp( -query => $query, -method_args => [ -num_alignments => 1 ]);               
    #         my $query_length = length($uniprots->{$annotn_uniprot_acc}->{seq});  
    #         my $hit = $results->next_hit();

    #         # check for results before going any further
    #         if ( !$results || !$hit) {
    #           next LINE;
    #         }
    #         my $hit_length = $hit->length();
    #         my $translation_stable_id = $hit->name();
    #         #my $evalue = $hit->hsp->evalue(); # DO WE WANT TO FILTER ON e-value too?           
    #         my $identity = $hit->hsp->percent_identity(); 

    #         my $max_length = ($hit_length, $query_length)[$hit_length < $query_length];
    #         my $min_length = ($hit_length, $query_length)[$hit_length > $query_length];

    #         my $query_length_covered = $min_length / $max_length ; # inforce 0 < ratio < 1

    #         print "IDENTITY : $identity ; query-target RATIO: $query_length_covered  \n";

    #         # when the dba is secondary we are more conservative. When it is the main dba (first one fetched/head of taxonomy branch) the values are reset to default

    #         if ($is_first_brch_dba eq 'true' ){ # When it is the first or only dba reset threshold values to default
    #           $current_identity_thr = $MIN_IDENTIY_THRESHOLD ;
    #           $current_query_length_thr = $QUERY_LENGTH_COVERED ;
    #         } else {
    #           $current_identity_thr = 100;
    #           $current_query_length_thr = 1;
    #         }

    #         if ( ($identity >= $current_identity_thr ) && ($query_length_covered >= $current_query_length_thr)  ) {
    #           print $output $result_line . ',1,' . $id_type . ',' . $branch_species . "," 
    #           . $translation_stable_id . "," . $identity . "\n" ;
    #         } 
    #       }
    #     }
    #   }
    #   $is_first_brch_dba = 'false';
    }


}

sub write_output {
    # my $self = shift;

    # #here, we flow out three parameters on branch 2: filename, size, and comp_size
    # $self->dataflow_output_id({'filename' => $self->param('filename'),
    #                            'orig_size' => $self->param('presize'),
    #                            'comp_size' => $self->param('postsize')},
    #                            2);

    # #note that there's also an event flown out on branch 1 by default, because
    # #we haven't overridden it.
}

sub post_cleanup {

}

# private method
# if the initial filename had a .gz
sub _create_post_operation_filename {
    # my ($pre_operation_filename, $gzip_flags) = @_;

    # if ($gzip_flags =~ /-d/) {
    #     return substr($pre_operation_filename, 0, -3);
    # } else {
    #     return $pre_operation_filename . ".gz";
    # }

}

sub find_translation { 
  my ( $self, $dba, $uniprot_acc, $locus_tag, $gene_name ) = @_;
  my $translation;
  my $identifier_type;
  my $dbentry_adaptor = $dba->get_adaptor("DBEntry");
  my @transcripts_ids = $dbentry_adaptor->list_transcript_ids_by_extids($uniprot_acc);
  my @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($uniprot_acc);

  if ( scalar(@gene_ids) == 0 ) {
    @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($locus_tag);
    $identifier_type = 'locus';
  }
  if ( scalar(@gene_ids) == 0 ) {
    @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($gene_name);
    $identifier_type = 'name';
  }
  
  my $translation_adaptor = $dba->get_adaptor("Translation");
  my $transcript_adaptor  = $dba->get_adaptor("Transcript");
  my $gene_adaptor        = $dba->get_adaptor("Gene");
  my $transcript;
  my $gene;

  if ( scalar(@transcripts_ids) >= 1 ) {
    my $transcript_id = $transcripts_ids[0];
    $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
    $translation = $translation_adaptor->fetch_by_Transcript($transcript);
    $identifier_type = 'accession';
  } elsif ( scalar(@gene_ids) >= 1 ) {
    $gene = $gene_adaptor->fetch_by_dbID( $gene_ids[0] );
    my @transcripts = @{ $transcript_adaptor->fetch_all_by_Gene($gene) };
    $transcript = $transcripts[0];
    $translation = $translation_adaptor->fetch_by_Transcript($transcript);
  }

  return $translation, $identifier_type;
} ## end sub find_translation

sub get_uniprot_seq {
  my ($self, $acc) = @_;

  my $search = $self->param('search');
  my $seq = $search->get('http://www.uniprot.org/uniprot/'.$acc.'.fasta');
  $seq =~ s/^>\S+\s+([^\n]+)\n//;
  my $des = $1;
  $seq =~ tr/\n//d;
  return {seq=>$seq,des=>$des};
}

sub get_uniparc_seq {
  my ($self, $acc) = @_;

  my $search =  $self->param('search');
  my $seq = $search->get('http://www.uniprot.org/uniparc/?query=' . $acc . '&columns=sequence&format=fasta');
  $seq =~ s/^>\S+\s+([^\n]+)\n//;
  my $des = $1;
  $seq =~ tr/\n//d;
  return {seq=>$seq,des=>$des};
}
1;
