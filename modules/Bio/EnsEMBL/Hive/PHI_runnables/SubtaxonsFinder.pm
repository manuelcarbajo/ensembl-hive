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
        'MAX_SUB_TAX_DBAS' => 6,
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
    my $annotn_tax_id = $self->param_required("species_tax_id");
    if (! $annotn_tax_id) {
         die 'annotn_tax_id is not defined'; # Will cause job to fail and leave a message in log_message
    } 

}


sub run {
    my $self = shift;
    my $annotn_tax_id = $self->param_required("species_tax_id");
    my $output_subtaxons =  $self->_get_subtaxons_dbas($annotn_tax_id);

    $self->param('subtaxons_dbas', $output_subtaxons);
}

=head2 _get_subtaxons_dbas
    
    Description: a private method that returns an arrayref of all dbadaptors of a supplied taxonomy_id ( = with its taxonomic children if any).
    In the eventuality of too many taxonomic children, the methood returns a maximum number of subtax dbas, limit defined by MAX_SUB_TAX_DBAS

=cut
sub _get_subtaxons_dbas {
  my ($self, $annotn_tax_id) = @_;
    my $lookup = $self->_get_lookup();
    my $branch_dbas = $lookup->get_all_by_taxon_branch($annotn_tax_id);
    my $nb_taxon_branch_dbas = scalar(@{$branch_dbas});
    my $MAX_SUB_TAX_DBAS = $self->param_required("MAX_SUB_TAX_DBAS");

    if ($nb_taxon_branch_dbas == 0) {
      print "\tNo dbs for this specie tax_id: $annotn_tax_id\n";
  } elsif ( $nb_taxon_branch_dbas > $MAX_SUB_TAX_DBAS ) {# if too many branch_dbas, limit the selection to MAX_SUB_TAX_DBAS
      print "\tToo many dbs for this specie (tax_id:$annotn_tax_id):  $nb_taxon_branch_dbas\n";
      $branch_dbas = $self->_limit_branch_dbas(@{$branch_dbas}); 
  } 

  return $branch_dbas;
}

=head2 _limit_branch_dbas
    
    Description: a private method that limits a supplied arrayref of dbas to a predefined number (MAX_SUB_TAX_DBAS)

=cut

sub _limit_branch_dbas {
  my ($self, @branch_dbas) = @_;
  my $new_branch_dbas;
  
  my $coredb_found = undef;
  my $nb_branches_picked = 0;
  my $i = 0;
  my $nb_taxon_branch_dbas = scalar(@branch_dbas);
  my $MAX_SUB_TAX_DBAS = $self->param_required("MAX_SUB_TAX_DBAS");

  # Go through all dbs and select just the coreDB
  while (!$coredb_found && $i < $nb_taxon_branch_dbas) {
    my $brch_dba = $branch_dbas[$i++]; 
    if ($brch_dba->dbc()->dbname() !~ /collection/) { 
      $coredb_found = 'true';
      push @{$new_branch_dbas}, $brch_dba;
      $nb_branches_picked ++;
     }
  }

  # Now add a limited number of collection DBs 
  $i=0;
  while ($nb_branches_picked < $MAX_SUB_TAX_DBAS  && $i < $nb_taxon_branch_dbas) {
    my $brch_dba = $branch_dbas[$i++]; 
    if ($brch_dba->dbc()->dbname() =~ /collection/) {
      push @{$new_branch_dbas}, $brch_dba;
      $nb_branches_picked ++;
     }
  }

  return $new_branch_dbas;
}

sub write_output {
    my ($self, @branch_dbas) = @_;
    my $fan_branch_code = $self->param('fan_branch_code');

    my $subtax_entries = $self->_build_output_hash();

    # "fan out" into fan_branch_code:
    $self->dataflow_output_id($subtax_entries, $fan_branch_code);
}



=head2 _build_output_hash 

    Description: a private method that returns a hash of parameters for all subjobs to fan. 
    Each job will have its own corresponding dba + the current job ($self) input fanned parameters from phiFileReader

=cut

sub _build_output_hash {
    my $self = shift;
    my @hashes = ();
    my $subtax_dbas = $self->param_required("subtaxons_dbas");
    # my $line_fields = [  'id_number', # 0
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
                         # 'core_db_host' # independent from $rows
                         # 'core_db_port' # independent from $rows


    #                     ---- NEW FIELDS ----
    #                    # '_species' # 15, specific species name of the subtaxon branch
    #                    # '_db_name' # 17 subtaxon db name
    #                   ];

    foreach my $dba (@$subtax_dbas) {
        my $job_param_hash = {};
        # for (my $i = 0 ; $i < scalar(@$line_fields); $i++)  {
        #     $job_param_hash->{ @$line_fields[ $i ] } = $self->param(@$line_fields[ $i ]);
        # }
        $job_param_hash->{ '_species' } = $dba->{ '_species' };
        $job_param_hash->{ '_dbname' } = $dba->dbc()->{'_dbname'};
        push @hashes, $job_param_hash;
    }

    return \@hashes;

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


1;
