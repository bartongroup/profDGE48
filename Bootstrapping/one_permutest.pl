#!/sw/bin/perl

=head1 NAME

B<one_permutest.pl>

=head1 DESCRIPTION

Permutation test for two samples. To be used by C<grid_launcher_DE.pl>.

=cut


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
my $stat = 'mean';
my $alpha = 0.5;
GetOptions (
  'infile|i=s' => \$infile,
  'outfile|o=s' => \$outfile,
  #'stat|s=s' => \$stat,
  #'alpha|a=f' => \$alpha,
  'nboot|n=i' => \$B,
  help => \$help,
  man => \$man, 
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;


die "Need -infile\n" unless defined $infile;
die "Need -outfile\n" unless defined $outfile;

print "Permutation test infile=$infile outfile=$outfile stat=$stat nboot=$B\n";

open IN, $infile or die "Cannot open batch file $infile\n";
open OUT, ">$outfile" or die "Huh? $!\n";
while(my $gene = <IN>)
{
  chomp $gene;
  my $s1 = <IN>; chomp $s1;
  my $s2 = <IN>; chomp $s2;

  my $x1 = pdl(split /\t/, $s1);
  my $x2 = pdl(split /\t/, $s2);
  my $m1 = mean($x1);
  my $m2 = mean($x2);
  
  my $p = PermutationTest($x1, $x2, $B, $stat, $alpha, 1);
  my $fc = log($m2 / $m1) / log(2);
  
  printf OUT "%s\t%.4g\t%.4g\t%.4g\t%.4g\n", $gene, $fc, $p, $m1, $m2;
} 
close IN;
close OUT;


=head1 OPTIONS

=over 5

=item B<-infile>=I<pathname>

Input file.

=item B<-outfile>=I<pathname>

Output file.

=item B<-nboot>=I<number>

Number of bootstraps. The default number is 5000000.

=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut

