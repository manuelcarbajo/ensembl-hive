=pod 

=head1 NAME

    Bio::EnsEMBL::Hive::PHI_runnables::xrefs_cleaner

=head1 SYNOPSIS

   
    
=head1 DESCRIPTION

    This is a generic RunnableDB module for cleaning all PHI-base associated xrefs from all the dbs in the registry for a particular division. 

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


package Bio::EnsEMBL::Hive::PHI_runnables::xrefs_cleaner;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Bio::EnsEMBL::LookUp::LocalLookUp;
use Scalar::Util qw(looks_like_number);

sub param_defaults {

    return {

        'core_db_host'       => 'mysql-eg-prod-2.ebi.ac.uk',# HOW TO MAKE THIS ONLY DEFAULT? MUST CHANGE IF PROVIDED WITH INPUT ARG!
        'core_db_port'       => 4289,
        'core_db_user'       => 'ensro',

        'tax_db_name'       => 'ensembl_compara_master', 
        'tax_db_host'       => 'mysql-eg-pan-prod.ebi.ac.uk',# HOW TO MAKE THIS ONLY DEFAULT? MUST CHANGE IF PROVIDED WITH INPUT ARG!
        'tax_db_port'       => 4276,
        'tax_db_user'       => 'ensro',
        'tax_dba_group'     => 'taxonomy',

        'MAX_SUB_TAX_DBAS' => 15,

        'fan_branch_code'   => 2,

    };
}

sub fetch_input {
    my $self = shift;


    my $core_db = $self->param_required('core_db_host');

    unless (defined $core_db) {
        die "core_db $core_db does not exist"; # Will cause job to fail and leave a message in log_message
    }

    my $core_db_port = $self->param_required('core_db_port');

    unless (defined $core_db_port) {
        die "core_db_port $core_db_port does not exist"; # Will cause job to fail and leave a message in log_message
    }
}


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

    

=cut

sub run {
    my $self = shift @_;

    my $specific_species= $self->param_required('_species');
    unless (defined $specific_species) {
        die "param specific_species does not exist"; # Will cause job to fail and leave a message in log_message
    }
    my $evidence = $self->param_required('_evidence');
    unless (defined $evidence) {
        die "param evidence does not exist"; # Will cause job to fail and leave a message in log_message
    }
    my $phi_entry = $self->param_required('phi_entry');
    unless (defined $phi_entry) {
        die "param phi_entry does not exist"; # Will cause job to fail and leave a message in log_message
    }
    print "phi_entry $phi_entry with $evidence for $specific_species\n ";

}

=head2 write_output

    Description : Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.
                  Here we rely on the dataflow mechanism to create jobs.

    param('fan_branch_code'): defines the branch where the fan of jobs is created (2 by default).

=cut

sub write_output {  # nothing to write out, but some dataflow to perform:
    my $self = shift @_;
}




1;
