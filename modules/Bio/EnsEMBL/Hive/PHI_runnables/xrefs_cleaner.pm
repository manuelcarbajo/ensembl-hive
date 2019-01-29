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

    };
}

sub fetch_input {
    my $self = shift;

    my $core_db_host = $self->param_required('core_db_host');    
    unless (defined $core_db_host) {
        die "core_db_host $core_db_host does not exist"; # Will cause job to fail and leave a message in log_message
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
     
    my $core_db_host = $self->param_required('core_db_host');  
    print "------- CLEAN_XREF ------------  $core_db_host   \n";
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
