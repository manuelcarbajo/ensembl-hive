=pod 

=head1 NAME

    Bio::EnsEMBL::Hive::PHI_runnables::protein_blaster

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


package Bio::EnsEMBL::Hive::PHI_runnables::protein_blaster;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use Bio::EnsEMBL::EGPipeline::Xref::BlastSearch;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Carp;



sub param_defaults {

    return {

        'core_db_host'       => 'mysql-eg-prod-2.ebi.ac.uk',# HOW TO MAKE THIS ONLY DEFAULT? MUST CHANGE IF PROVIDED WITH INPUT ARG!
        'core_db_port'       => 4289,
        'core_db_user'       => 'ensro',

        'fan_branch_code'   => 2,
    };
}

sub fetch_input {
    my $self = shift;
    
    my $phi_entry = $self->param_required('phi_entry');    
    unless (defined $phi_entry) {
        die "phi_entry $phi_entry does not exist"; # Will cause job to fail and leave a message in log_message
    }

    my $core_db = $self->param_required('core_db_host');

    unless (defined $core_db) {
        die "core_db $core_db does not exist"; # Will cause job to fail and leave a message in log_message
    }

    my $core_db_port = $self->param_required('core_db_port');

    unless (defined $core_db_port) {
        die "core_db_port $core_db_port does not exist"; # Will cause job to fail and leave a message in log_message
    }

    my $annotn_uniprot_acc = $self->param_required('uniprot_acc');
    if (! $annotn_uniprot_acc || _remove_spaces($annotn_uniprot_acc) eq '' ) {
        die 'uniprot accesion is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $branch_species = $self->param_required('_species');
    if (! $branch_species || _remove_spaces($branch_species) eq '' ) {
        die 'Branch species is not defined'; # Will cause job to fail and leave a message in log_message
    }

    my $db_name = $self->param_required("_dbname");
    if (! $db_name) {
         die "dbname is not defined"; # Will cause job to fail and leave a message in log_message
    }
}


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;


    my $annotn_uniprot_acc = $self->param_required('uniprot_acc');
    my $branch_species = $self->param_required('_species');
    my $core_db = $self->param_required('core_db_host');
    my $core_db_port = $self->param_required('core_db_port');
    my $db_name = $self->param_required("_dbname");

    my $brch_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new( -host    => $core_db ,
                                                        -port    => $core_db_port ,
                                                        -user    => 'ensro',
                                                        -dbname  => $db_name,
                                                        -species => $branch_species,
                                                        -multispecies_db => 1
                                                       );

    my ($uniprots, $id_type) = $self->get_sequence($annotn_uniprot_acc);

    if ($uniprots->{$annotn_uniprot_acc}->{seq} ne '') { 
        print "\t-Got uniprot/uniparc accession.- \n";


    #  Create BLAST DB for this specie if it doesn't already exists or if it needs to be updated

    my $first_char = substr($branch_species, 0, 1); # to avoid too many DBs in a folder we subdivide by qalphabetical order
    my $blast_db_dir = $self->param('blast_db_directory') . "$first_char" ;
    my $blast_database_file = $blast_db_dir . "/" . $branch_species . ".pep.all.fa";
    my $checksum_file = $blast_db_dir . "/" . $branch_species . ".cks";
    print "BLAST search: creating peptide file if it doesnt exist already in $blast_database_file:\n";


    my $pept_file_update_needed = $self->is_update_needed($brch_dba, $branch_species);

    
    if ( ! -f $blast_database_file || $pept_file_update_needed) { # Returns undef if the file doesn't exist, or if the file exists but it's not plain file (e.g. if it's a directory).
      open (my $db_file, ">", $blast_database_file) or croak "Could not open ".$blast_database_file;;
      print "update_needed:$pept_file_update_needed, getting all genes for $branch_species \n";

      my $translation_adaptor=$brch_dba->get_TranslationAdaptor();
      my @translations = @{$translation_adaptor->fetch_all()}; 

      foreach my $translation (@translations) {
        my $stable_id = $translation->stable_id();
        my $sequence = $translation->seq();
        print $db_file ">" . $stable_id . "\n";
        print $db_file $sequence . "\n";
      }
      close $db_file;

      my $checksum = $self->param('checksum');

      open(my $cksfile, '>', $checksum_file);
      print $cksfile $checksum;
      close $cksfile;
    
    } else {
      print "Update not needed, using previously existing DB.\n";
    }

    #     # Create a BlastPlus factory
    #     my $fac;
       
    #     eval {
    #       $fac = Bio::Tools::Run::StandAloneBlastPlus->new( 
    #         -db_name => $branch_species,
    #         -db_dir => $blast_db_dir ,
    #         -db_data => $blast_database_file,
    #         -create => 1
    #       );
    #       $fac->make_db();
    #     };
    #     $@ =~ /EXCEPTION|WARNING/ and my $e = $@;
    #     if ( defined $e ) {
    #       print "\t******************* WARNING: **************************\n";
    #       print "BlastPlus factory could not be created for $annotn_uniprot_acc: This entry will be skipped\n";
    #       next LINE;
    #     }

    #     my $query   = Bio::Seq->new( -display_id => $annotn_uniprot_acc, -seq =>  $uniprots->{$annotn_uniprot_acc}->{seq} );
    #     my $results = $fac->blastp( -query => $query, -method_args => [ -num_alignments => 1 ]);               
    #     my $query_length = length($uniprots->{$annotn_uniprot_acc}->{seq});  
    #     my $hit = $results->next_hit();

    #     # check for results before going any further
    #     if ( !$results || !$hit) {
    #       next LINE;
    #     }
    #     my $hit_length = $hit->length();
    #     my $translation_stable_id = $hit->name();
    #     #my $evalue = $hit->hsp->evalue(); # DO WE WANT TO FILTER ON e-value too?           
    #     my $identity = $hit->hsp->percent_identity(); 

    #     my $max_length = ($hit_length, $query_length)[$hit_length < $query_length];
    #     my $min_length = ($hit_length, $query_length)[$hit_length > $query_length];

    #     my $query_length_covered = $min_length / $max_length ; # inforce 0 < ratio < 1

    #     print "IDENTITY : $identity ; query-target RATIO: $query_length_covered  \n";

    #     # when the dba is secondary we are more conservative. When it is the main dba (first one fetched/head of taxonomy branch) the values are reset to default

    #     if ($is_first_brch_dba eq 'true' ){ # When it is the first or only dba reset threshold values to default
    #       $current_identity_thr = $MIN_IDENTIY_THRESHOLD ;
    #       $current_query_length_thr = $QUERY_LENGTH_COVERED ;
    #     } else {
    #       $current_identity_thr = 100;
    #       $current_query_length_thr = 1;
    #     }

    #     if ( ($identity >= $current_identity_thr ) && ($query_length_covered >= $current_query_length_thr)  ) {
    #       print $output $result_line . ',1,' . $id_type . ',' . $branch_species . "," 
    #       . $translation_stable_id . "," . $identity . "\n" ;
    #     } 
    }    

}

=head2 is_update_needed

    Description : private method to determine if a peptide file for a given species need to be created, updated or none.

=cut

sub is_update_needed   {   
    my ($self, $brch_dba, $branch_species) = @_;

    my $blast_db_dir = $self->param('blast_db_directory') . substr($branch_species, 0, 1) ;
    my $checksum_file = $blast_db_dir . "/" . $branch_species . ".cks";

    my $dbc=$brch_dba->dbc();
    my $temp_table_name = "temp_prots_$branch_species";
    my $sth = $dbc->prepare( "CREATE TEMPORARY TABLE $temp_table_name select count(*) 
                                from translation tl join transcript tc using(transcript_id) 
                                join seq_region sr using(seq_region_id) 
                                join coord_system cs using(coord_system_id) 
                                join meta m using(species_id) 
                                where meta_key='species.production_name' and meta_value like '$branch_species';");
    $sth->execute() or
     die "Error creating temporary table for $branch_species: perhaps the DB doesn't have a meta table or the table already exists?\n" .
       "$DBI::err .... $DBI::errstr\n";

    $sth = $dbc->prepare( "checksum TABLE $temp_table_name ;");
    $sth->execute() or
     die "Error obtaining checksum for temporary table $temp_table_name\n" .
       "$DBI::err .... $DBI::errstr\n";

    my $table_name;
    my $checksum;
    $sth->bind_columns(\$table_name,\$checksum);
    $sth->fetch();
    print "Checksum for $branch_species = $checksum\n";
    $self->param('checksum',$checksum);
    
    # Move the delete statement to a bulk loop at the end of db_load_hasher
    $sth = $dbc->prepare( "DROP TABLE $temp_table_name;");
    $sth->execute();
    
    my $pept_file_update_needed = undef;

    if ( ! -f $checksum_file ) {
      $pept_file_update_needed = 1;
    } else {
        open my $file, '<', $checksum_file; 
        my $stored_checksum = <$file>; 
        close $file;
        if ($stored_checksum ne $checksum ) {
            print "$checksum is different than expected checksum $stored_checksum. Blast peptide fasta file update needed\n";
            $pept_file_update_needed = 1;
        }
    }
    
    return $pept_file_update_needed;
}

=head2 write_output

    Description : Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.
                  Here we rely on the dataflow mechanism to create jobs.

    param('fan_branch_code'): defines the branch where the fan of jobs is created (2 by default).

=cut

sub write_output {  # nothing to write out, but some dataflow to perform:
    my $self = shift @_;
    #TODO
    print "This vaues should now go to the accumulator\n";
}

=head2 get_sequence 

    Description: a private method that returns the protein sequence corresponding to a given accession.
    Insure that uniprot/uniparc queries only one protein sequence and no more. 

=cut

sub get_sequence {
  my ($self, $annotn_uniprot_acc) = @_;

  my $id_type;
  my $uniprots = {};
  my $up;

  eval {
    $up = $self->get_uniprot_seq($annotn_uniprot_acc);
    $uniprots->{$annotn_uniprot_acc} = $up;
    $id_type = 'uniprot';
  };

  if($@ || $uniprots->{$annotn_uniprot_acc}->{seq} eq '') {
    $up = $self->get_uniparc_seq($annotn_uniprot_acc);
    $uniprots->{$annotn_uniprot_acc} = $up;
    $id_type = 'uniparc';
  }

  # Insure that uniprot/uniparc queries only one protein sequence and no more. 
  my $count_seqs = $up->{seq}   =~ tr/>//; 
  if ($count_seqs > 0) {
    my $annotn_phi_entry = $self->param_required('phi_entry');
    die  " *-*-*- ERROR : skipping $annotn_phi_entry  *-*-* Uniprot/Uniparc query returns more than one sequence\n";
  }

  return ($uniprots, $id_type);
}

sub get_uniprot_seq {
  my ($self, $acc) = @_;
  
  my $search = Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();
  my $seq = $search->get('http://www.uniprot.org/uniprot/'.$acc.'.fasta');
  $seq =~ s/^>\S+\s+([^\n]+)\n//;
  my $des = $1;
  $seq =~ tr/\n//d;
  return {seq=>$seq,des=>$des};
}

=head2 get_uniprot_seq 

    Description: a private method that returns a protein sequence (or more!) from a accession query to UNIPROT

=cut

sub get_uniprot_seq {
  my ($self, $acc) = @_;
  
  my $search = Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();
  my $seq = $search->get('http://www.uniprot.org/uniprot/'.$acc.'.fasta');
  $seq =~ s/^>\S+\s+([^\n]+)\n//;
  my $des = $1;
  $seq =~ tr/\n//d;
  return {seq=>$seq,des=>$des};
}

=head2 get_uniparc_seq 

    Description: a private method that returns a protein sequence (or more!) from a accession query to UNIPARC

=cut
sub get_uniparc_seq {
  my ($self, $acc) = @_;
  
  my $search =  Bio::EnsEMBL::EGPipeline::Xref::BlastSearch->new();
  my $seq = $search->get('http://www.uniprot.org/uniparc/?query=' . $acc . '&columns=sequence&format=fasta');
  $seq =~ s/^>\S+\s+([^\n]+)\n//;
  my $des = $1;
  $seq =~ tr/\n//d;
  return {seq=>$seq,des=>$des};
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
