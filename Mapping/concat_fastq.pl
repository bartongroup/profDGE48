#!/usr/bin/perl

=head1 NAME

concat_fastq.pl - script to concatenate technical replicates into biological replicates

=cut

use strict;
use warnings;

use Getopt::Long qw(:config auto_version);
use Pod::Usage;
use File::Basename;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $inpath;
my $outpath;
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
our $VERSION = '0.1';

GetOptions (
   'src-path=s'=> \$inpath,
   'dest-path=s' => \$outpath,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid source path.') unless ($inpath);
pod2usage(-msg => 'Please supply a valid destination path.') unless ($outpath);

## check that the destination has 'Snf2' and 'WT' directories
foreach my $d (qw/Snf2 WT/) {
   die "ERROR - destination directory doesn't have a '$d/' subdirectory\n" unless ( -d "$outpath/$d" );
}

## read source directory with the assumed structure as below
# src-path/
#   1/
#     120903_0219_D0PT7ACXX_1_SA-PE-096.sanfastq.gz
#     120903_0219_D0PT7ACXX_2_SA-PE-096.sanfastq.gz
#     etc.
#   2/
#     120903_0219_D0PT7ACXX_1_SA-PE-021.sanfastq.gz
#     120903_0219_D0PT7ACXX_2_SA-PE-021.sanfastq.gz
#     etc.
#   ..
#   96/
#     120903_0219_D0PT7ACXX_1_SA-PE-085.sanfastq.gz
#     120903_0219_D0PT7ACXX_2_SA-PE-085.sanfastq.gz
#     etc.
#
## subdirs 1-48 are assumed to be 'Snf2' and 49-96 are 'WT'
print "Reading $inpath...\n" if $VERBOSE;
opendir(my $top, $inpath) or die "ERROR - unable to open directory '$inpath': $!\nDied";
while (my $entry = readdir($top)) {
   next if ($entry =~ /^\./);
   if (-d "$inpath/$entry") {
      next unless ($entry =~ /^\d+$/); # skip any non-numeric directory names
      #print "Found $entry/\n";
      opendir(my $sub, "$inpath/$entry") or die "ERROR - unable to open directory '$inpath/$entry': $!\nDied";
      my $mid;
      my @files;
      
      # open subdirs, get MID number and store filenames
      while (my $subentry = readdir($sub)) {
         next if ($subentry =~ /^\./);
         # exctract MID from filename
         if ($subentry =~ /(\d{3}).sanfastq.gz$/) {
            $mid = sprintf("%02d",$1);
         } else {
            die "ERROR - unable to extract MID from '$subentry' filename\nDied"
         }
         push @files, "$inpath/$entry/$subentry";
      }
      closedir($sub);
      
      # assume there must be seven files per dir
      die "ERROR - wrong number of files found in $entry/ (".scalar @files."). Expected 7.\n" unless (scalar @files == 7);
      
      # do 'cat'ing...
      my $CMD = '';
      if ($entry < 49) {  # first 48 are Snf2
         my $rep = sprintf("%02d",$entry);
         $CMD = "zcat @files > $outpath/Snf2/Snf2_rep${rep}_MID${mid}_allLanes.fastq";
      } else {            # and second 48 are WT
         my $rep = sprintf("%02d",$entry-48);
         $CMD = "zcat @files > $outpath/WT/WT_rep${rep}_MID${mid}_allLanes.fastq";
      }
      die "ERROR - unable to create command entry '$entry'\nDied" unless (defined($CMD));
      print "Concatenating files for sample $entry...\n" if $VERBOSE;
      print "CMD: $CMD\n" if $DEBUG;
      system("$CMD") == 0 or die "ERROR - system() failed for '$CMD'\nDied";
   }
}
closedir($top);
print "Finished\n";

=head1 SYNOPSIS

concat_fastq.pl --src-path <path> --dest-path <path> [--version] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

This script takes the multilane, technical replicate data (as produced by GenePool) and combines it into biological replicate fastq files.

It assumes the source directory has 96 subdirs (numbered 1-96) and that the first 48 are 'Snf2' samples and the second 48 are 'WT' samples. It also assumes that in each of the number subdirs there are seven technical replicate files, which need combining into one, biological, replicate.

Output files are put in the destination directory under 'Snf2' and 'WT' subdirs. Filenames are derived from the biological replicate number, sample and MID. e.g.:

   WT_rep06_MID02_allLanes.fastq 

=head1 OPTIONS

=over 5

=item B<--src-path>

Source directory with all replicate data. Must have numeric subdirs with fastq files.

=item B<--dest-path>

Destination directory for saving concatenated fastq files

=item B<--version>

Report version info and exit

=item B<--verbose|--no-verbose>

Toggle verbosity. [default:none]

=item B<--debug|--no-debug>

Toggle debugging output. [default:none]

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back

=head1 AUTHOR

Chris Cole <christian@cole.name>

=head1 COPYRIGHT

Copyright 2012, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut