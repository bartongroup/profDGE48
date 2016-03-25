#!/sw/bin/perl 

=head1 NAME

B<one_bootstrap_ET93.pl>

=head1 DESCRIPTION

Bootstrap test for two samples, to be run by grid launcher. Following Efron & Tibshirani (1993), "An introduction to the bootstrap", p.220-224.

To be used by C<grid_launcher_DE.pl> and it is not supposed to be ran on its own.

=cut


#
# Tested against t-test for normal samples, 30/10/2012
# Modified and tested again, 1/04/2013
#

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;

use CompBio::Statistics;

$| = 1;

my ($infile, $outfile);
my ($help, $man);
my $B = 5000000;
GetOptions (
  'infile|i=s' => \$infile,
  'outfile|o=s' => \$outfile,
  'nboot|n=i' => \$B,
  help => \$help,
  man => \$man, 
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;

die "Need -infile\n" unless defined $infile;
die "Need -outfile\n" unless defined $outfile;

print "Bootstrap with:\n";
print "infile $infile\n";
print "outfile $outfile\n";
print "nboot $B\n";

open IN, $infile or die "Cannot open batch file $infile\n";
open OUT, ">$outfile" or die "Huh? $!\n";
while(my $gene = <IN>)
{
  chomp $gene;
  my $s1 = <IN>; chomp $s1;
  my $s2 = <IN>; chomp $s2;

  my $x1 = pdl(split /\t/, $s1);
  my $x2 = pdl(split /\t/, $s2);
  
  print $gene;
  my ($p, $m1, $m2) = BootstrapTest($x1, $x2, $B, 1);
  print "\n";
  my $fc = log($m2 / $m1) / log(2);
  
  printf OUT "%s\t%.4g\t%.4g\t%.4g\t%.4g\n", $gene, $fc, $p, $m1, $m2;
} 
close IN;
close OUT;





=head1 OPTIONS

=over 5

=item B<-infile>=I<pathname>

Name of the input file.

=item B<-outfile>=I<pathname>

Name of the output file.

=item B<-nboot>=I<number>

Number of bootstraps. The default value is 5000000.

=back

=cut
