#!/sw/bin/perl -w

=head1 NAME

B<count_fastq_reads_lanes.pl> - count reads in each condition, replicate and lane. Creates a file ''all_readcount_lane.dat'' with counts. Requires directory whith subdirs, one per biological sample, containing lean FASTQ files.

=cut


use strict;
use warnings;

use Getopt::Long;

use PDL;
use PDL::NiceSlice;

use GRNASeq;

$| = 1;


my $fastqdir;
GetOptions(
  'dir=s' => \$fastqdir
);

#my $fastqdir = "$topdir/raw/raw_reads";
my $fileformat = '';


open F, ">all_readcount_lane.dat" or die;
for my $cond (@conds)
{
  for my $rep (1 .. 48)
  {
    for my $lane (1 .. 7)
    {
      my $file = GetFile($cond, $rep, $lane);
      print "$cond $rep $lane ";
      my $n = CountReads($file);
      print F "$cond $rep $lane $file $n\n";
      print "$n\n";
    }
  }
}

########################################################

sub GetFile
{
  my ($cond, $rep, $lane) = @_;
  
  my $sample = ($cond eq 'Snf2') ? $rep : $rep + 48; 
  
  my $dir = "$fastqdir/$sample";
  die unless -d $dir;
  opendir D, $dir or die;
  my $file;
  while(my $f = readdir D)
    {$file = "$dir/$f" if $f =~ /ACXX_$lane/}
  return $file 
}

sub CountReads
{
  my $file = shift;
  my @s = `gunzip -c $file | wc`;
  my ($lines) = split " ", $s[0];
  return $lines / 4
}


=head1 SYNOPSIS

  count_fastq_reads_lanes.pl -dir=/dir/with/fastq/files
      
=head1 OPTIONS

=over 4

=item B<-dir>=I<path>

Directory with FASTQ files. File names should be in a specific format, so this very simple script really works only with our experiment.


=back
