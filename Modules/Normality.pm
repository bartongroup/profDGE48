package CompBio::Normality;

=head1 NAME

CompBio::Normality - normality test

=head1 SYNOPSIS

  use PDL;
  use CompBio::Normality;
  
  my $x = grandom(10);
  my $P = NormalityTest($x);
  
=head1 DESCRIPTION
  
Test for sample normality. Based on D.Agostino, Belanger & D'Agostino, "A suggestion for using powerful and informative tests of normality", The American Statistician, 44, 316.
    
=head1 FUNCTIONS

=over 4

=cut

require Exporter;

use strict;
use warnings;

use PDL;
use PDL::NiceSlice;
use Math::CDF qw(:all);
use Carp;

our $VERSION = 1.0;

our @ISA = ("Exporter");
our @EXPORT = qw(
  Skewness
  NormalityTest
  SkewnessTest
  KurtosisTest
);

################################################################

=item Skewness

  $sk = Skewness($d);

  Unbiased estimator of sample skewness.

=cut


sub Skewness
{
  my $d = shift;
  
  my ($m, $s) = stats($d);
  my $n = nelem($d);
  carp "n must be greater than 3" if $n <= 3;
  
  my $m2 = sum(($d - $m)**2) / $n;
  my $m3 = sum(($d - $m)**3) / $n;

  return 0 if $m2 == 0;
  my $g1 = $m3 / $m2**1.5;
  my $cor = sqrt($n*($n - 1)) / ($n - 2);  # unbias
  
  return $cor * $g1
}

=item NormalityTest

  $p = NormalityTest($d);
  
Test that the underlying distribution is normal. $d is a sample (one-dimensional piddle). Returns the corresponding p-value. When $p is small, the underlying distribution is likely to be not normal.

This is the omnibus test based on skewness and curtosis.
  
=cut

sub NormalityTest
{
  my ($d) = @_;
  
  return 1 if nelem($d) <= 3;
  my $Z1 = SkewnessTest($d);
  my $Z2 = KurtosisTest($d);
  my $K2 = $Z1*$Z1 + $Z2*$Z2;
  
  my $p = 1 - pchisq($K2, 2);
}

sub SkewnessTest
{
  my ($d) = @_;
  
  my ($m, $s) = stats($d);
  my $n = nelem($d);
  carp "n must be greater than 3" if $n <= 3;
  
  my $m2 = sum(($d - $m)**2) / $n;
  my $m3 = sum(($d - $m)**3) / $n;

  my $sb1 = $m3 / $m2**1.5;
  my $Y = $sb1 * sqrt(($n + 1) * ($n + 3) / (6 * ($n - 2)));
  my $b2 = 3 * ($n*$n + 27*$n - 70)*($n + 1)*($n + 3) / (($n - 2)*($n + 5)*($n + 7)*($n + 9));
  my $W = sqrt(-1 + sqrt(2 * ($b2 - 1)));
  $W += 1e-6 if $W == 1;                 # can I do this???
  my $delta = 1 / sqrt(log($W));
  my $alpha = sqrt(2 / ($W**2 - 1));
  my $Z = $delta * log($Y/$alpha + sqrt(($Y/$alpha)**2 + 1));
}

sub KurtosisTest
{
  my ($d) = @_;
  
  my ($m, $s) = stats($d);
  my $n = nelem($d);
  carp "n must be greater than 3" if $n <= 3;
  
  my $m2 = sum(($d - $m)**2) / $n;
  my $m4 = sum(($d - $m)**4) / $n;
  my $b2 = $m4 / $m2**2;
  my $Eb2 = 3*($n - 1) / ($n + 1);
  my $Vb2 = 24 * $n * ($n - 2)*($n - 3) / (($n + 1)**2 * ($n + 3)*($n + 5));
  my $x = ($b2 - $Eb2) / sqrt($Vb2);
  my $sb12 = 6*($n*$n - 5*$n + 2) / (($n + 7)*($n + 9)) * sqrt(6*($n + 3)*($n + 5) / ($n*($n - 2)*($n - 3)));
  my $A = 6 + 8/$sb12 * (2/$sb12 + sqrt(1 + 4/$sb12**2));
  my $Z = ((1 - 2/(9*$A)) - ((1 - 2/$A) / (1 + $x*sqrt(2/($A - 4))))**(1/3))/sqrt(2 / (9*$A));
}

1;
