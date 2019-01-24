=pod 

=head1 NAME

Bio::EnsEMBL::Hive::PHI_runnables::SubtaxonsFinder

=head1 SYNOPSIS

    standaloneJob.pl Bio::EnsEMBL::Hive::PHI_runnables::SubtaxonsFinder -input_id_list '[{"a"=>1},{"a"=>2}]' -debug 1

=head1 DESCRIPTION

    This runnable [de]compresses a particular file, and makes note of the size
    of that file before and after [de]compression

    This operation could also be done by chaining together steps using the default
    SystemCmd Runable; here they are all combined into a single Runnable
    as a demonstration.

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


package Bio::EnsEMBL::Hive::PHI_runnables::TranslationFinder;

use strict;
use warnings;
use Bio::EnsEMBL::EGPipeline::Xref::BlastSearch;
use Bio::EnsEMBL::LookUp::LocalLookUp;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::EGPipeline::Xref::BlastSearch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::LookUp::LocalLookUp;
use Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor;


use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Time::HiRes qw( time );

use Scalar::Util qw(reftype);

sub param_defaults {

        return {
        'ready_branch_code' => 3,
        'blast_branch_code' => 4,
    };
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
    # working, etc. 
  
    my $annotn_tax_id = $self->param_required("species_tax_id");
    if (! $annotn_tax_id) {
         die 'annotn_tax_id is not defined'; # Will cause job to fail and leave a message in log_message
    } 

    my $specific_species_name = $self->param_required("_species");
    if (! $specific_species_name) {
         die "specific_species_name ('_species')is not defined"; # Will cause job to fail and leave a message in log_message
    }
    
    my $dbname = $self->param_required("_dbname");
    if (! $dbname) {
         die "dbname is not defined"; # Will cause job to fail and leave a message in log_message
    }

}


sub run {
    my $self = shift;
    my $specific_species_name = $self->param_required("_species");
    my $translation;
    my $id_type;
    my $brch_dba;
    ($translation, $id_type, $brch_dba) = $self->_find_translation();
    
    $self->param('translation', $translation);
    $self->param('id_type', $id_type);
    $self->param('brch_dba', $brch_dba);
  
}

sub write_output {
    my $self = shift;

    my $blast_branch_code = $self->param('blast_branch_code');
    my $ready_branch_code = $self->param('ready_branch_code');
    my $translation = $self->param('translation');

    my $id_type = $self->param('id_type');
    my $brch_dba = $self->param('brch_dba');

    my $translation_found_params ;
    my $translation_unknown_params ;

    if($translation) {
      print "Translation_found:" . $translation->stable_id . ":\n";
      $translation_found_params = $self->_build_output_hash($translation, $id_type, $brch_dba);
      print Dumper($translation_found_params);
      $self->dataflow_output_id($translation_found_params, $ready_branch_code);
      
    } else {
      print "Translation_unknw:$translation:\n";
      $translation_unknown_params = $self->_build_output_hash($translation, $id_type, $brch_dba);
      print Dumper($translation_unknown_params);
      $self->dataflow_output_id($translation_unknown_params, $blast_branch_code);
    }
    
}



=head2 _get_lookup
    
    Description: a private method that loads the registry and the lookup taxonomy DB.

=cut

sub _get_lookup {
    my $self = shift @_;
    Bio::EnsEMBL::Registry->load_registry_from_db( 
                         -USER => 'ensro',
                         -HOST => 'mysql-eg-prod-2.ebi.ac.uk',
                         -PORT => 4239);
      
    
    my $lookup =      
                Bio::EnsEMBL::LookUp::LocalLookUp->new( -SKIP_CONTIGS => 1,
                                                        -NO_CACHE     => 1 ); # Used to check all leafs sub_specie/strains under a taxonomy_id (specie)

    my $tax_dba_details = undef;


    $tax_dba_details->{-GROUP}  = 'taxonomy';
    $tax_dba_details->{-DBNAME} = 'ncbi_taxonomy';
    $tax_dba_details->{-HOST}   = 'mysql-eg-staging-2.ebi.ac.uk';
    $tax_dba_details->{-PORT}   = 4275;
    $tax_dba_details->{-USER}   = 'ensro';

    my $tax_adaptor = Bio::EnsEMBL::DBSQL::TaxonomyNodeAdaptor->new(
                            Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{$tax_dba_details} ) 
                      );
    
    $lookup->taxonomy_adaptor($tax_adaptor);

    return ($lookup);
}

=head2 _find_translation
    
    Description: a private method that tries to find a translation in the Ensemble DB by either protein_accession, locus or by gene_type.

=cut

sub _find_translation { 
  my $self = shift;
  
  my $specific_species_name = $self->param_required("_species");
  my $specific_db_name = $self->param_required("_dbname");
  print "-----------------------------  $specific_species_name in $specific_db_name  ---------------------------------\n";
  my $core_db = $self->param_required("core_db_host");
  my $core_db_port = $self->param_required("core_db_port");
  my $lookup = $self->_get_lookup();
  
  #build dbadaptor from specific_db_name and specific_species_name. This is messy but nI didn't find another way of passing the dba along the analysis.
  #If you find a better way of doing please improve this
  my $dba =
    Bio::EnsEMBL::DBSQL::DBAdaptor->new( -host    => $core_db ,
                                         -port    => $core_db_port ,
                                         -user    => 'ensro',
                                         -dbname  => $specific_db_name,
                                         -species => $specific_species_name,
                                         -multispecies_db => 1
                          );
  my $dbc=$dba->dbc();
  my $sth = $dbc->prepare( "SELECT DISTINCT species_id FROM $specific_db_name.meta " .
        "WHERE meta_key='species.production_name' AND meta_value LIKE '$specific_species_name'");
  $sth->execute() or
    die "Error querying for species_id: perhaps the DB doesn't have a meta table?\n" .
      "$DBI::err .... $DBI::errstr\n";

  my $species_id;
  $sth->bind_columns(\$species_id);
  $sth->fetch;
  
  $dba->species_id($species_id);
  # end get $dba... Oufff...!!
  
  my $uniprot_acc = $self->param_required('uniprot_acc');
  my $locus_tag = $self->param_required('locus');
  my $gene_name = $self->param_required('gene_name');
  print "-----------------------------  accession $uniprot_acc /locus $locus_tag /gene_name $gene_name ---------------------------------\n";
  # if ($specific_species_name eq 'fusarium_graminearum_gca_000599445') {
  #   print "dba:\n";
  #   print Dumper ($dba );
  # }
  my $translation = '';
  my $identifier_type = '';

  my $dbentry_adaptor = $dba->get_adaptor("DBEntry");
  my @transcripts_ids = $dbentry_adaptor->list_transcript_ids_by_extids($uniprot_acc); 

  print " nb transcript ids by uniprot:" . scalar(@transcripts_ids) . "\n";
  my @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($uniprot_acc); # List of gene_id by an external identifier accession that is
                                                                          # linked to  any of the genes transcripts, translations or the gene itself
  print " nb gene_ids by uniprot:" . scalar(@gene_ids) . "\n";                                                                          

  if ( scalar(@gene_ids) == 0 ) {
    @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($locus_tag);   # List of gene_id by an external identifier locus that is
                                                                          # linked to  any of the genes transcripts, translations or the gene itself
    $identifier_type = 'locus';                                               

    print " nb gene_ids by locus:" . scalar(@gene_ids) . "\n";
  }
  if ( scalar(@gene_ids) == 0 ) {
    @gene_ids = $dbentry_adaptor->list_gene_ids_by_extids($gene_name);   # List of gene_id by an external identifier gene name that is
                                                                          # linked to  any of the genes transcripts, translations or the gene itself
    $identifier_type = 'name';
    print " nb gene_ids by name:" . scalar(@gene_ids) . "\n";
  }

  my $translation_adaptor = $dba->get_adaptor("Translation");
  my $transcript_adaptor  = $dba->get_adaptor("Transcript");
  my $gene_adaptor        = $dba->get_adaptor("Gene");
  my $transcript;
  my $gene;
  
  if ( scalar(@transcripts_ids) >= 1 ) {
    
    my $transcript_id = $transcripts_ids[0];
    print "CONTROL++++++ ++a+ uniprot_acc: $uniprot_acc:, locus_tag: $locus_tag:, gene_name: $gene_name:, identifier_type: accession; transcript_id: $transcript_id\n";
    $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
    $translation = $translation_adaptor->fetch_by_Transcript($transcript);
    $identifier_type = 'accession';
  } elsif ( scalar(@gene_ids) >= 1 ) {
      print "CONTROL++++++ ++b+gene_ids[0]:" . $gene_ids[0] . ":, uniprot_acc: $uniprot_acc:, locus_tag: $locus_tag:, gene_name: $gene_name:, identifier_type: $identifier_type;\n";
      $gene = $gene_adaptor->fetch_by_dbID( $gene_ids[0] );
      my @transcripts = @{ $transcript_adaptor->fetch_all_by_Gene($gene) };

      $transcript = $transcripts[0];
      $translation = $translation_adaptor->fetch_by_Transcript($transcript);
  } else {
      print "NO TRANSCRIPTS NO GENES FOR THIS:: uniprot_acc: $uniprot_acc:, locus_tag: $locus_tag:, gene_name: $gene_name:, identifier_type: $identifier_type;\n";
  }
  print "CONTROL++++++ END find_translation+++:\n";
  return $translation, $identifier_type, $dba;
} 

=head2 _build_output_hash 

    Description: a private method that returns a hash of parameters to add to the input plus() 

=cut

sub _build_output_hash {
    my ($self, $translation, $id_type, $dba) = @_;

        # line_fields =    'id_number', # 0
    #                      'phi_entry', # 1
    #                      'uniprot_acc', # 2
    #                      'gene_name', # 3
    #                      'locus', # 4
    #                      'origin', # 5
    #                      'species_tax_id', # 6
    #                      'subtaxon_id', # 7
    #                      'pathogen_name', # 8
    #                      'host_tax_id', # 9
    #                      'host_name', # 10
    #                      'phenotype', # 11
    #                      'experiment_condition', # 12
    #                      'litterature_id', # 13
    #                      'DOI' # 14
    #                      '_species' # 15, specific species name of the subtaxon branch
    #                      '_db_name' # 16 subtaxon db name
    #
    #                       -------- NEW FIELDS -----------
    #                      
    #                      '_evidence' # 17 is this a direct_match or a blast_match
    #                      '_id_type' # 18 How the match was found (if DIRECT_MATCH locus_tag, gene name or accession), or Protein DB used to get the sequence blast_match
    #                      '_translation_id' # 19 translation stable _id
    #                      '_percent_identity' # 20 percentage identity of the match
    #

    my $job_param_hash = {};
    my $evidence = 'BLAST_MATCH';

    if ($translation) {
      print "\t\t--- Found gene on ensembl: ";
      my $gene_adaptor        = $dba->get_adaptor("Gene");
      my $gene = $gene_adaptor->fetch_by_translation_stable_id($translation->stable_id);
          
      print "using $id_type\n\t\t" . $gene->stable_id . "\t (should have a translation):" . $translation->stable_id . ":\n";
      
      $evidence = 'DIRECT_MATCH';
      $job_param_hash->{ '_id_type' } = $id_type;
      $job_param_hash->{ '_translation_id' } = $translation->stable_id;
      $job_param_hash->{ '_percent_identity' } = 100;

    } else {
        print "\t\tDidn't find direct match on Ensembl \t id_type :$id_type\n";

        $job_param_hash->{ '_id_type' } = undef;
        $job_param_hash->{ '_translation_id' } = undef;
        $job_param_hash->{ '_percent_identity' } = undef;
      }

    $self->param('evidence', $evidence);
    $job_param_hash->{ '_evidence' } = $evidence ;

    return $job_param_hash;
  
}

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
