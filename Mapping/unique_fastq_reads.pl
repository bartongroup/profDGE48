#!/usr/bin/perl

=head1 NAME

unique_fastq_reads.pl - create a Fastq file with only unique examples of reads

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $fastq;
my $out = 'unique.fastq';
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
my $version = 0.1;

GetOptions (
   'in=s'      => \$fastq,
   'out=s'     => \$out,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid filename.') unless ($fastq && -s $fastq);

print "Running version: $version\n\n";
print "Reading '$fastq'...\n" if $VERBOSE;
my $header = '';
my $seq = '';
my $length;
my %uniqReads;
my $n = 0;
my $z = new IO::Uncompress::Gunzip $fastq or die " ERROR - gunzip failed on '$fastq': $GunzipError\n";
while (my $line = $z->getline()) {
   chomp($line);
   if ($. % 4 == 1) {  # read header
      $header = $line;
      ++$n;
   }
   if ($. % 4 == 2) {  # sequence
      $seq = $line;
   }
   if ($. % 4 == 0) {  # qulaity string
      $uniqReads{$seq} = "$header|$line";
      $seq = '';
      $header = '';
   }
}
$z->close();
printf "Read $n reads of which %d are unique\n", scalar keys %uniqReads if $VERBOSE;

print "Writing output to '$out'...\n" if $VERBOSE;
open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
foreach my $s (keys %uniqReads) {
   my ($head,$qual) = split(/\|/, $uniqReads{$s});
   print $OUT "$head\n$s\n+\n$qual\n";
}
close($OUT);
print "Done!\n" if $VERBOSE;

=head1 SYNOPSIS

unique_fastq_reads.pl --in <file> [--out <file>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

Takes a Fastq file and writes out only unique examples of reads therein. Keeps header and quality information for the first instance of each read.

=head1 OPTIONS

=over 5

=item B<--in>

Input fastq file (can be gzipped).

=item B<--out>

Output filename. [default: unique.fastq]

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