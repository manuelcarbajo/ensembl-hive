=pod 

=head1 NAME

    Bio::EnsEMBL::Hive::PHI_runnables::phi_accumulator_writer

=head1 SYNOPSIS

   
    
=head1 DESCRIPTION

    This is a module for validating all the phi_base candidate entries and grouping them into a super_hash in the accumulator. 

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


package Bio::EnsEMBL::Hive::PHI_runnables::phi_accumulator_writer;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;


sub param_defaults {

    return {

        'tax_db_name'       => 'ensembl_compara_master', 
        'tax_db_host'       => 'mysql-eg-pan-prod.ebi.ac.uk',# HOW TO MAKE THIS ONLY DEFAULT? MUST CHANGE IF PROVIDED WITH INPUT ARG!
        'tax_db_port'       => 4276,
        'tax_db_user'       => 'ensro',
        'tax_dba_group'     => 'taxonomy',

    };
}

sub fetch_input {
    my $self = shift;

    my $_species = $self->param_required('_species');    
    unless (defined $_species) {
        die "_species $_species does not exist"; # Will cause job to fail and leave a message in log_message
    }

    my $_dbname = $self->param_required('_dbname');      
    unless (defined $_dbname) {
        die "_dbname $_dbname does not exist for $_species"; # Will cause job to fail and leave a message in log_message
    }

    my $phibase_id = $self->param_required('phi_entry');
    if (! $phibase_id || _remove_spaces($phibase_id) eq '' ) {
         die 'phibase_id is not defined'; # Will cause job to fail and leave a message in log_message
    } 
    
    my $translation_stable_id = $self->param_required('_translation_id');
    if (! $translation_stable_id || _remove_spaces($translation_stable_id) eq '' ) {
         die " translation_stable_id is not defined for $phibase_id _dbname $_dbname"; # Will cause job to fail and leave a message in log_message
    }

    my $acc = $self->param_required('uniprot_acc');
    if (! $acc || _remove_spaces($acc) eq '' ) {
         die 'uniprot accesion is not defined'; # Will cause job to fail and leave a message in log_message
    } 

    my $host_tax_id = $self->param_required('host_tax_id');
    if (! $host_tax_id || _remove_spaces($host_tax_id) eq '' ) {
         die ' host_tax_id is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $host_name = $self->param_required('host_name');
    if (! $host_name || _remove_spaces($host_name) eq '' ) {
         die ' host_name is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $phenotype_name = $self->param_required('phenotype');
    if (! $phenotype_name || _remove_spaces($phenotype_name) eq '' ) {
         die ' phenotype_name is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $evidence = $self->param_required('_evidence');
    if (! $evidence || _remove_spaces($evidence) eq '' ) {
         die " evidence is not defined for $phibase_id _dbname $_dbname"; # Will cause job to fail and leave a message in log_message
    }
       
    print "PHI ACCUMULATOR HASH translation: $translation_stable_id for $phibase_id // $_dbname  ::: $_species :: $evidence \n";
}


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
     

}

=head2 write_output

    Description : Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.
                  Here we rely on the dataflow mechanism to create jobs.

    param('fan_branch_code'): defines the branch where the fan of jobs is created (2 by default).

=cut

sub write_output {  # nothing to write out, but some dataflow to perform:
    my $self = shift @_;
}

=head2 _remove_spaces 

    Description: a private method that returns a string without leading or trailing whitespaces

=cut
sub _remove_spaces {
  my $string = shift;
  $string =~ s/^\s*//;
  $string =~ s/\s*$//;
  return $string;
}


1;
