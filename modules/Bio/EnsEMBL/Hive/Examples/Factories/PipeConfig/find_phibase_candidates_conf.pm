=pod

=head1 NAME

  Bio::EnsEMBL::Hive::Examples::Factories::PipeConfig::find_phi_base_candidates_conf

=head1 SYNOPSIS


    init_pipeline.pl Bio::EnsEMBL::Hive::Examples::Factories::PipeConfig::find_phi_base_candidates_conf  -pipeline_url $EHIVE_URL -hive_force_init 1

    seed_pipeline.pl -url $EHIVE_URL -logic_name inputfile -input_id "{ 'inputfile' => '/homes/mcarbajo/phytopath_dbi/rothamstead_small_version.csv' }"

    runWorker.pl -url 'mysql://ensrw:scr1b3d1@mysql-eg-devel-1.ebi.ac.uk:4126/mcarbajo_ehive_phibase_db' -debug 1


=head1 DESCRIPTION

    This is an example pipeline put together from five basic building blocks:

    Analysis_1: JobFactory.pm is used to turn the list of files in a given directory into jobs

        these jobs are sent down the branch #2 into the second analysis

    Analysis_2: JobFactory.pm is used to run a wc command to determine the size of a file
                (and format the output with sed), then capture the command's object, putting
                the file size into a parameter for later use.

    Analysis_3: SystemCmd.pm is used to run these compression/decompression jobs in parallel.

    Analysis_4: JobFactory.pm is used to run a wc command to determine the size of a file
                (and format the output with sed), then capture the command's object, putting
                the file size into a parameter for later use.

    Analysis_5: SystemCmd.pm is used to run the notify-send command, displaying a message on the screen.

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

package Bio::EnsEMBL::Hive::Examples::Factories::PipeConfig::find_phibase_candidates_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf; # Need to use as well as use base() to turn on INPUT_PLUS()


#Defines which parameters are required from the user command's line
sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
     %{$self->SUPER::pipeline_wide_parameters},
    'inputfile'     => $self->o('inputfile'),
  };
}

#defines default values for some of the parameters
sub default_options {
  my ($self) = @_;
  return {
    %{$self->SUPER::default_options(@_)},

    #number_of_taxa_to_consider => '10',# just an example
  };
}

=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.
                  Here it defines five analyses:

                    * 'find_files'          generates a list of files whose names match the pattern #only_files#
                                            Each job of this analysis will dataflow (create jobs) via branch #2 into 'compress_a_file' analysis.

                    * 'pre_compress_size'   determines the size of a file, and captures the size in a parameter to pass through the pipeline

                    * 'compress_a_file'     actually performs the (un)zipping of the files in parallel

                    * 'post_compress_size'   determines the size of a file, and captures the size in a parameter to pass through the pipeline

                    * 'notify'              displays a message on the screen using the notify-send command

=cut

sub pipeline_analyses {
    my ($self) = @_;
    return [
            {
                -logic_name => 'inputfile',
                -module     => 'Bio::EnsEMBL::Hive::PHI_runnables::phiFileReader',
                -parameters => {
                    'inputfile' => '~/phytopath_dbi/rothamstead_small_version.csv',
                    'delimiter' => ',',
                    'column_names' => 1,
                    'output_ids' => '#output_ids#',
                },
                -flow_into => {
                # Create a fan of jobs, using INPUT_PLUS() to propagate all of
                # this analysis' parameters down branch 2
                2 => { 'find_subtaxons' => INPUT_PLUS() },
                },

            },
            {
                -logic_name    => 'find_subtaxons',
                -module        => 'Bio::EnsEMBL::Hive::PHI_runnables::SubtaxonsFinder', 
            },
           ];
   
}

1;