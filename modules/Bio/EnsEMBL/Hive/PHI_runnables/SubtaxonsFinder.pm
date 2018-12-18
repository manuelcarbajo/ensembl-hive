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


package Bio::EnsEMBL::Hive::PHI_runnables::SubtaxonsFinder;

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
  
        'delimiter'         => ',',
        'lookup'            => undef,

        'tax_db_name'       => 'ensembl_compara_master', 
        'tax_db_host'       => 'mysql-eg-pan-prod.ebi.ac.uk',# HOW TO MAKE THIS ONLY DEFAULT? MUST CHANGE IF PROVIDED WITH INPUT ARG!
        'tax_db_port'       => 4276,
        'tax_db_user'       => 'ensro',
        'tax_dba_group'     => 'taxonomy',
        'tax_db_species' => undef,
        'tax_db_pass' => undef,
        'tax_db_driver' => undef,
        'tax_db_species_id' => 1,
        'tax_db_multispecies_db' => 0,

        'MAX_SUB_TAX_DBAS' => 15,

        'fan_branch_code'   => 2,
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

    my $annotn_tax_id = $self->param_required("_6");
    if (! $annotn_tax_id) {
         die 'annotn_tax_id is not defined'; # Will cause job to fail and leave a message in log_message
    }

    
    # # use param_required for this, since we want the job to fail if there's no output_ids.
    # my $output_ids = $self->param_required('input_id');


    # if (! $output_ids) {
    #     die 'output_ids is not defined'; # Will cause job to fail and leave a message in log_message
    # }

}


sub run {
    my $self = shift;
    
    
    my $annotn_tax_id = $self->param_required("_6");

    my $search = Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();
    $self->param('search', $search);
    my $lookup = _get_lookup();
    my $branch_dbas = $lookup->get_all_by_taxon_branch($annotn_tax_id);
    my $nb_taxon_branch_dbas = scalar(@{$branch_dbas});

    for my $brch_dba (@{$branch_dbas}) {

      my $mc = $brch_dba->get_MetaContainer();
      my $branch_species = $mc->single_value_by_key('species.production_name');
      my $division =  $mc->single_value_by_key('species.division');

      print " ----  branch_species: $branch_species \n";
    }


}

sub write_output {
    # my $self = shift;

my $fan_branch_code         = $self->param('fan_branch_code');
}

sub post_cleanup {

}

=head2 _get_lookup
    
    Description: a private method that loads the registry and the lookup taxonomy DB.

=cut

sub _get_lookup {
    my $self = shift @_;
    Bio::EnsEMBL::Registry->load_registry_from_db( 
                         -USER => 'ensro',
                         -HOST => 'mysql-eg-prod-1.ebi.ac.uk',
                         -PORT => 4238);
      
    
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
