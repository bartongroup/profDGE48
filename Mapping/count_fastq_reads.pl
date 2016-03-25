#!/sw/bin/perl -w

=head1 NAME

B<count_fastq_reads.pl> - count reads in all fastq files and store results in output files <cond>_readcount.dat  

=cut


use strict;
use warnings;

use PDL;
use PDL::NiceSlice;
use Getopt::Long;

use GRNASeq;

$| = 1;


my $fastqdir;
GetOptions(
  'dir=s' => \$fastqdir
);

my $fileformat = '';



for my $cond (@conds)
{
  my %files = GetFileList($cond);
  open F, ">${cond}_readcount.dat" or die;
  for my $rep (sort {$a <=> $b} keys %files)
  {
    print "$cond$rep $files{$rep}";
    my $n = CountReads($files{$rep});
    print F "$files{$rep} $rep $n\n";
    print "   $n\n";
  }
}

sub GetFileList
{
  my ($cond) = @_;
  
  my $dir = "$fastqdir/$cond";
  die unless -d $dir;
  opendir D, $dir or die;
  my %files = ();
  while(my $f = readdir D)
    {$files{int($1)} = "$dir/$f" if $f =~ /${cond}_rep(\d+)_MID(\d+)/}
  return %files 
}

sub CountReads
{
  my $file = shift;
  my @s = `gunzip -c $file | wc`;
  my ($lines) = split " ", $s[0];
  return $lines / 4
}


=head1 SYNOPSIS

  count_fastq_reads.pl -dir=/dir/with/fastq/files
      
=head1 OPTIONS

=over 4

=item B<-dir>=I<path>

Directory with FASTQ files. File names should be in a specific format, so this very simple script really works only with our experiment.


=back