#!/usr/bin/perl

=head1 NAME

summary_stats.pl - generate summary statistics of alignment infomation

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use threads;
use threads::shared;
use Thread::Queue;

my $path;
my $samtools = '/sw/samtools-0.1.18/samtools';
my $out = 'out.tsv';
my $numThreads = 1;
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
my $VERSION = '0.2';

GetOptions (
   'path=s'    => \$path,
   'samtools=s' => \$samtools,
   'out=s'     => \$out,
   'threads=i' => \$numThreads,
   'version!'  =>  sub { print "$0 version: $VERSION\n"; exit },
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid path that is readable.') unless ($path && -d $path && -r $path);

## output version info
print "\n$0 version: $VERSION\n";
printf "$samtools version: %s\n", getVersion($samtools);

## get tophat paths
print "Searching $path for tophat directories...\n" if $VERBOSE;
opendir(my $dh, $path) or die "ERROR - unable to open directory '$path': $!\nDied";
my $queue = Thread::Queue->new();
while (my $entry = readdir($dh)) {
   next if ($entry =~ /^\./);        # skip hidden/special files
   next unless (-d "$path/$entry");  # skip anything that isn't a directory
   
   $queue->enqueue("$path/$entry");  # add dirs to thread queue
}
closedir($dh);
printf "Found %d directories to read\n", $queue->pending() if $VERBOSE;
my $ndir = $queue->pending();

## check don't have too many threads
if ($numThreads > $queue->pending()) {
   printf "DEBUG Reducing number of threads to: %d\n", $queue->pending() if $DEBUG;
   $numThreads = $queue->pending();
}

## get read counts for mapped and unmapped reads
print "Getting read counts from tophat runs...\n"  if $VERBOSE;
my %mapped :shared;
my %unmapped :shared;
my %failed :shared;
my @threads;
for (my $i = 0; $i < $numThreads; ++$i) {
   $threads[$i] = threads->create(\&getReadCounts, \%mapped, \%unmapped, \%failed);
}
$_->join foreach (@threads);
printf "Found %d mapped files and %d unmapped files\n", scalar keys %mapped, scalar keys %unmapped if $VERBOSE;
warn "Warning - mismatch in number of tophat directories and the number of mapped files\n", if (scalar keys %mapped != $ndir);
warn "Warning - mismatch in number of tophat directories and the number of unmapped files\n", if (scalar keys %unmapped != $ndir);

open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
print $OUT "Directory\tMappedReads\tUnmappedReads\tFailedReads\n";
foreach my $p (sort keys %mapped) {
   my $base = basename($p);
   print $OUT "$base\t$mapped{$p}\t$unmapped{$p}\t$failed{$p}\n";
}

# get the version string for a given program name
sub getVersion {
   my $prog = shift;
   
   my $version = '';
   if ($prog =~ /samtools/) {
      my $out = `$prog 2>&1`;
      if ($out =~ /Version: ([\.\d]+)/) {
         $version = $1;
      }
   } elsif ($prog =~ /htseq-count/) {
      my $out = `$prog 2>&1`;
      if ($out =~ /version (.*)\./) {
         $version = $1;
      }
   } else {
      die "ERROR - program '$prog' not known. Can't extract version number from it\nFix the code, man!\nDied";
   }
   
   if ($version) {
      return($version)
   } else {
      warn "Warning - unable to get version info for '$prog'\n";
      return(0)
   }
}

## get read counts for tophat produced BAM files
sub getReadCounts {
   my $mapped = shift;
   my $unmapped = shift;
   my $failed = shift;
   
   while(my $path = $queue->dequeue()) {
      foreach my $file (qw/accepted_hits.bam unmapped.bam/) {
         print "Getting stats for $path/$file...\n";
         if (!-s "$path/$file") {
            warn "Warning - file '$path/$file' doesn't exist or is empty. Skipping\n";
            next;
         }
         my $out = `$samtools flagstat $path/$file` or die "ERROR - unable to run samtools: $!\nDied";
         if ($file eq 'accepted_hits.bam') {
            if ($out =~ /(\d+).*mapped/) {
               $mapped{$path} = $1;
               #print "Mapped: $1\n";
               #printf "No entries %d\n", scalar keys %$mapped;
            } else {
               die "ERROR - unable to retrieve mapped read counts for '$path/$file'\nDied";
            }
            
         } else {
            if ($out =~ /(\d+) \+ (\d+) in total/) {
               $unmapped{$path} = $1;
               $failed{$path} = $2;
               #print "Unmapped: $1\n";
            }
         }
      }
      return(0) if ($queue->pending() == 0);
   }
}

=head1 SYNOPSIS

summary_stats.pl --path <path> [--samtools <path>] [--threads <num>] [--out <file>] [--version] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

This script generates simple statistics on the outcome of Tophat alignments.

It expects that the given I<--path> includes subdirectories generated by tophat which include I<accepted_hits.bam> and I<unmapped.bam> files in them.

Given that assumption the script with then run 'samtools flagstat' on the BAM files to determine the number of mapped, unmapped and failed reads for each run. These are then reported in the output file.

=head1 OPTIONS

=over 5

=item B<--path>

Path to tophat output dirs.

=item B<--samtools>

Path to samtools executeable [default: /sw/samtools-0.1.18/samtools]

=item B<--threads>

Number of threads to use [default: 1]

=item B<--out>

Output filename. [default: out.tsv]

=item B<--version>

Output version info and exit.

=item B<--verbose|--no-verbose>

Toggle verbosity. [default: verbose]

=item B<--debug|--no-debug>

Toggle debugging output. [default: no-debug]

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