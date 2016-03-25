#!/sw/bin/perl -w

=head1 NAME

B<make_bamlinks.pl> - Creates symbolic links to .bam files in subdirectories.

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

my $topdir;
my ($help, $man);
GetOptions(
  'topdir=s' => \$topdir,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
die "Need -topdir\n" unless defined $topdir;
die "Cannot find $topdir\n" unless -d $topdir;

chdir $topdir;
opendir D, "." or die;
my @lst = readdir D;
for my $d (@lst)
{
  next unless -d $d;
  next if $d eq '.' || $d eq '..';
  my $bam = "$d/accepted_hits.bam";
  my $link = "$d.bam";
  unlink $link;
  die "BAM file not found in $d\n" unless -e $bam;
  `ln -s $bam $link\n`;
}


=head1 SYNOPSIS

  make_bamlinks.pl -topdir=/cluster/gjb_lab/cdr/GRNAseq/mapping/genome/WT

=head1 DESCRIPTION

C<run_alignements.pl> and C<run_biological_alignements.pl> scripts create individual directories for each replicate under one top (condition) directory. Each directory contain C<accepted_hits.bam> file. However, further analysis requires that all bam files are in one directory. This script creates symbolic links from the top directory to individual bam files.

=head1 OPTIONS

=over 5

=item B<-topdir>

Top directory containing bam subdirectories. Typically, this contains one condition.

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back

=cut