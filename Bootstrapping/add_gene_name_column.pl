#!/usr/bin/perl

=head1 NAME

add_gene_name_column.pl - given a tab delimited file with Ensembl IDs in the first col, add a gene name.

=cut

use strict;
use warnings;

use Getopt::Long qw(VersionMessage :config auto_version);
use Pod::Usage;
use File::Basename;

#use lib '/sw/opt/ensembl-api/bioperl-live';
#use lib '/sw/opt/ensembl-api/68/ensembl/modules';
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;

my $file;
my $species = '';
my $feature = 'Gene';
my $desc = 0;
my $coords = 0;
my $out = 'ensembl_annotated.csv';
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
our $VERSION = "1.4";

$|=1;

GetOptions (
   'in=s'      => \$file,
   'species=s' => \$species,
   'feature-type=s' => \$feature,
   'desc!'     => \$desc,
   'coords!'   => \$coords,
   'out=s'     => \$out,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
   'version'   => sub { VersionMessage( -msg => sprintf("ensembl API %s\n", software_version()) ) }
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid filename.') unless ($file && -e $file);
pod2usage(-msg => 'Please supply a species name.') unless ($species);

$desc = 1 if ($coords); # enable descriptions if co-ordinates are set.

die "ERROR - invalid feature-type '$feature'. Use 'Gene' or 'Transcript'." unless ($feature =~ /^(gene|transcript)$/i); 

$species = lc($species);
$species = 's_cerevisiae' if ($species eq 'yeast');

printf "NOTE: using Ensembl API version %s\n", software_version();
my $registry = connectEnsemblRegistry($species);
my $adaptor = $registry->get_adaptor($species, 'Core', 'Gene');
die "ERROR - failed to get adaptor for '$species'. Check spelling and that it's a valid Ensembl species. Or check that you're using the correct API.\n" unless (defined($adaptor));
warn "Warning - API version check has failed. You probably need to update your local install.\n" unless ($registry->version_check($registry->get_DBAdaptor($species, 'core')));

print "Annotating entries...\n" if $VERBOSE;
open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
open(my $IN, "<", $file) or die "ERROR - unable to open '$file': ${!}\nDied";
my $c = 0;
while(<$IN>) {
   chomp;
   s/\"//g; # remove quotes;
   my @F = split(/\t/, $_);
   my $id = shift @F;
   if ($id =~ /^GeneID/) { # print header line
      if ($desc) {
         printf $OUT "$id\t%s\tGeneName\tDescription",join("\t",@F);
         if ($coords) {
            print $OUT "\tChromosome\tStart\tEnd\tStrand\n"
         } else {
            print $OUT "\n";
         }
      } else {
         printf $OUT "$id\t%s\tGeneName\n",join("\t",@F);
      }
      next;
   }
   ++$c;
   print "Fetching gene '$id' from Ensembl...\n" if $DEBUG;
   my $g;
   if ($feature =~ /transcript/i) {
      $g = $adaptor->fetch_by_transcript_stable_id($id)
   } else {
      $g = $adaptor->fetch_by_stable_id($id); # fetch it from ensembl
   }
   my $name;
   my $gDesc;
   my $chr = '-';
   my $start = '-';
   my $end = '-';
   my $strand = '-';
   
   if (!defined($g)) {  # check gene exists
      warn "Warning - gene '$id' not found\n";
      $name = 'unknown';
      $gDesc = 'none';
   } else {
      $name = $g->external_name(); # get the name
      $name = 'unknown' unless ($name);
      if ($desc) {
         $gDesc = $g->description();
         $gDesc = 'none' unless ($gDesc);
         $gDesc =~ s/\[Source.*\]//;
      }
      if ($coords) {
         $chr = $g->slice->seq_region_name();
         $start = $g->start();
         $end = $g->end();
         $strand = $g->strand();
      }
      
   }
   if ($desc) {
      printf $OUT "$id\t%s\t$name",join("\t",@F);
      print $OUT "\t$gDesc";  # print $gDesc separately as can sometimes have unexpected chars (e.g. %) which confuse printf()
      if ($coords) {
         print $OUT "\t$chr\t$start\t$end\t$strand\n";
      } else {
         print $OUT "\n";
      }
   } else {
      printf $OUT "$id\t%s\t$name\n",join("\t",@F);
   }
   
}
close($IN);
if ($c == 0) {
   warn "Warning - no gene names found. Nothing changed\n";
   unlink($out);
   exit;
}
#rename($out, $fasta) or die "Failed to rename '$out' to '$fasta'\n";
print "\nDone!\n";

## this is required in order to pick the correct
## connection parameters to the Ensembl API as 
## species from the Ensembl Genomes projects differ from the main API
sub connectEnsemblRegistry {
   my $species = shift;
   
   my $registry = 'Bio::EnsEMBL::Registry';
   my %main;
   $main{$_}++ foreach (qw/chicken human mouse s_cerevisiae/);
   
   if ($main{$species}) {  # this is for the main API species
      print "Connect to main Ensembl API...\n" if $VERBOSE;
      
      $registry->load_registry_from_db(
          -host => 'ensembldb.ensembl.org',
          -user => 'anonymous'
      );
      
   } else {  # this is for the Ensemble Genomes API species
      print "Connecting to Ensembl Genomes API...\n" if $VERBOSE;
      
      $registry->load_registry_from_db(
          -host => 'mysql.ebi.ac.uk',
          -user => 'anonymous',
          -port => 4157,
      );
      
   }
   return($registry);
}



=head1 SYNOPSIS

add_gene_name_column.pl --in <file> --species <name> [--feature-type <string>] [--desc|--no-desc] [--coords|--no-coords] [--out <file>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help] [--version]

=head1 DESCRIPTION

Given a tab-delimited list with Ensembl gene IDs in the first column, add gene names (and optionally descriptions and coordinates) to the list.

The addtional information is placed as the last set of columns in the data, leaving the original data at the same columns.

If the file has a header line, the first column must be called 'GeneID' otherwise the formatting will be wrong.

=head1 OPTIONS

=over 5

=item B<--in>

Input tab-delimited file. 1st column must be an ensembl ID.

=item B<--species>

Species name.

=item B<--feature-type>

Specify whether the Ensembl IDs correspond to genes or transcripts. [default: Genes]

=item B<--desc|--no-desc>

Toggle whether to add gene description as well as gene name. [default: off]

=item B<--coords|--no-coords>

Toggle whether to include genomic coordinates (assumes --desc is set). [default: off]

=item B<--out>

Output filename. [default: ensembl_annotated.csv]

=item B<--verbose|--no-verbose>

Toggle verbosity. [default: on]

=item B<--debug|--no-debug>

Toggle debugging output. [default: off]

=item B<--version>

Only print out the version information before exiting.

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back

=head1 AUTHOR

Chris Cole <christian@cole.name>

=cut