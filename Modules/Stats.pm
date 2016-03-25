package Stats;
require Exporter;
require PDL;
require PDL::NiceSlice;
require Math::CDF;
require Carp;
require Gamma;
require Tools;

=head1 NAME

CompBio::Statistics - a few statistical subroutines

=head1 SYNOPSIS

  use CompBio::Statistics;
  
  my ($a, $b, $SEa, $SEb) = LinFit($x, $y); # linear fit
  my $r = corr($x, $y); # Pearson's correlation coefficient
  my $r = ucorr($x, $y); # unbiased estimator of Pearson's correlation coefficient
  my $g = genGauss($mean, $sigma);   # generate number with Gaussian distribution
  my $p = genPoisson($mean);         # generate number with Poissonian distribution
  my ($D, $Prob) = KStwo($data1, $data2);  # Two-sample KS test
  my ($p, $t, $dof) = tTest ($d1, $d2);    # t-test for samples with identical variances
  my ($p, $t, $dof) = tuTest ($d1, $d2);   # t-test for samples with different variances
  my ($median, $ml, $mu) = MedianCI($x, $gamma) # median confidence interval
  my $pc = BenjaminiHochberg($p); # Benjamini-Hochberg correction
    
=head1 FUNCTIONS

=over 4

=cut


use strict;
use warnings;

#use Math::Trig;
use PDL;
use PDL::NiceSlice;
use Math::CDF qw(:all);
use Carp;
use Gamma;
use Tools;

our $VERSION = 2.3;

our @ISA = ("Exporter");
our @EXPORT = qw(
 Gammaln
 Betai
 Betacf
 genGauss
 genPoisson
 MeanStddev
 LinFit
 LinFit0
 wLinFit
 mean
 corr
 corrb
 ucorr
 rank
 Spearman
 KSone
 KStwo
 KSone2D
 quaduni1
 tTest
 tuTest
 logtTest
 logntTest
 LogrankTest
 OneSampletTest
 ShrinktTest
 OneSampleShrinktTest
 ShrinkVar
 ShrinkVarUnitTest
 MannWhitney
 PermutationTest
 BootstrapTest
 TrimmedMean
 EnrichmentTest
 UnitTestEnrichmentTest
 InitLogfact
 hypergeom
 fast_hypergeom
 binomial
 MedianCI
 PropCI
 CohensKappa
 QKS
 Disagreement
 $Betacf_warning
);

our $Betacf_warning;
my $PI = 4 * atan(1);

################################################################

=item Gammaln

  my $g = Gammaln($z);

  Log of gamma function, ln(Gamma(z)), z > 0 (from Numerical Recipes). Not very
  accurate. Please use CompBio::SpecFun::Gamma instead.

=cut


sub Gammaln
{
  my $z = shift;

  my @cof = (76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5);
  
  my $x = $z - 1;
  my $tmp = $x + 5.5;
  $tmp = ($x + 0.5)*log($tmp) - $tmp;
  my $ser = 1;
  for my $c (@cof)
  {
    $x += 1;
    $ser += $c / $x;
  }
  return $tmp + log(2.5066282746310005 * $ser);
}

#########################################################

sub Betai
# Incomplete beta function I_x(a,b)
{
  my ($a, $b, $x) = @_;
  my $bt;

  #print "a = $a, b = $b, x = $x in Betai\n";
  die "Bad x = $x in routine Betai\n" if ($x < 0.0 || $x > 1.0 );
  if ($x == 0.0 || $x == 1.0 )
    {$bt = 0.0}
  else
    {$bt= exp(Gammaln($a+$b) - Gammaln($a) - Gammaln($b) + $a * log($x) + $b * log(1.0 - $x))}

  if($x < ($a + 1.0)/($a + $b + 2.0))
    {return $bt * Betacf($a, $b, $x) / $a}
  else
    {return 1.0- $bt * Betacf($b, $a, 1.0-$x) / $b}
}

#########################################################

sub Betacf
{
  my ($a,  $b,  $x) = @_;

  #undef $Betacf_warning;

  my $MAXIT = 10000;
  my $EPS = 3.0e-7;
  my $FPMIN = 1.0e-80;

  my $qab = $a + $b;
  my $qap = $a + 1.0;
  my $qam = $a - 1.0;
  my $c = 1.0;
 
  my $d = 1.0 - $qab * $x/$qap;
  $d = $FPMIN if abs($d) < $FPMIN;
 
  $d = 1.0/$d;
  my $h = $d;
  my $del = 2;
  for my $m (1 .. $MAXIT)
  {
    my $m2 = 2*$m;
    my $aa = $m*($b-$m)*$x/(($qam + $m2) * ($a+$m2));
    $d = 1.0 + $aa*$d;
    $d = $FPMIN if abs($d) < $FPMIN;
    $c = 1.0 + $aa/$c;
    $c = $FPMIN if abs($c) < $FPMIN;
    $d = 1.0/$d;
    $h *= $d*$c;
    $aa = -($a+$m) * ($qab + $m)*$x/(($a + $m2) * ($qap + $m2));
    $d = 1.0 + $aa * $d;
    $d = $FPMIN if abs($d) < $FPMIN;
    $c = 1.0 + $aa/$c;
    $c = $FPMIN if abs($c) < $FPMIN;
    $d = 1.0/$d;
    $del = $d*$c;
    $h *= $del;
    last if abs($del-1.0) > $EPS;
  }
  if (abs($del-1.0) <= $EPS)
  {
    #warn "a = $a b = $b x = $x del = $del h = $h\n";
    #warn "a or b too big, or MAXIT too small in Betacf\n";
    #warn "No convergence in Betafc\n";
    #warn "Don't trust your results!\n";
    $Betacf_warning = 1;
    $h = 1 if $h > 1;
    $h = 0 if $h < 0;
  }
  return $h;
}

################################################################

=item genGauss

  my $g = genGauss($mean, $sigma);

Generates a number with Gaussian distribution of a given mean and standard deviation, using Box-Muller method.

=cut

sub genGauss
{
  my ($mean, $sigma) = @_;
  
  return $mean if $sigma == 0;
  my ($x1, $x2, $w, $y1);
  do{{
    $x1 = 2 * rand() - 1;
    $x2 = 2 * rand() - 1;
    $w = $x1 * $x1 + $x2 * $x2;
  }} while($w >= 1);
  $w = sqrt((-2 * log($w))/$w);
  $y1 = $x1 * $w;
  
  $y1 = $y1 * $sigma + $mean;
  return($y1);
}

################################################################

sub MeanStddev
#
#  ($mean, $sigma) = MeanStddev(\@lc);
#
# Finds mean and standard deviation for a given array of data,
# using a canning little recursive algorithm allowing for one loop and
# avoiding potential dangers of rounding up on subtraction and
# producing negative variance.
{
  my $dat = shift;
  
  my $m = $dat->[0];
  my $s = 0;
  my $p = 1;
  for my $x (@$dat)
  {
    if($p > 1)
    {
      $m = (($p - 1)*$m + $x)/$p;
      $s = (($p - 2)*$s + $p/($p - 1) * ($x - $m)**2)/($p - 1);
    }
    $p++
  }
  $s = sqrt($s);
  return($m, $s);
}

######################################################

sub LinFit
#
# y = ax + b
#
{
  my ($x, $y) = @_;
  
  my $sel = which(($x->isgood()) & ($y->isgood()));
  $x = $x($sel);
  $y = $y($sel);
  
  return (0, 0) if nelem($x) < 3 || nelem($y) < 3 || nelem($x) != nelem($y);
  my $n = nelem($x);
  my $Mx = sum($x) / $n;
  my $My = sum($y) / $n;
  my $Sxx = sum(($x - $Mx)**2);
  my $Syy = sum(($y - $My)**2);
  my $Sxy = sum(($x - $Mx) * ($y - $My));
  
  return (0, 0) if $Sxx == 0;
  
  my $a = $Sxy / $Sxx;
  my $b = $My - $a * $Mx;

  my $SDR = sqrt(($Syy - $a * $Sxy) / ($n - 2));
  my $SEa = $SDR / sqrt($Sxx);  # standard error
  my $SEb = $SDR * sqrt(1/$n + $Mx**2 / $Sxx);

  return ($a, $b, $SEa, $SEb, $Mx, $Sxx, $SDR);
}

sub LinFit_
#
# y = ax + b
#
{
  my ($x, $y) = @_;
  
  return (0, 0) if nelem($x) < 3 || nelem($y) < 3 || nelem($x) != nelem($y);
  my $n = nelem($x);
  my $Sx = sum($x);
  my $Sy = sum($y);
  my $Sxx = sum($x * $x);
  my $Syy = sum($y * $y);
  my $Sxy = sum($x * $y);
  my $SSxx = $Sxx - $Sx * $Sx / $n;
  my $SSyy = $Syy - $Sy * $Sy / $n;
  my $SSxy = $Sxy - $Sx * $Sy / $n;
  return (0, 0) if $SSxx == 0;
  
  my $a = $SSxy / $SSxx;
  my $b = ($Sy - $a * $Sx) / $n;

  my $s = sqrt(($SSyy - $a * $SSxy) / ($n - 2));
  my $SEb = $s * sqrt(1 / $n + ($Sx/$n)**2 / $SSxx);  # standard error
  my $SEa = $s / sqrt($SSxx);

  return ($a, $b, $SEa, $SEb);
}

######################################################

sub LinFit0
#
# y = ax
#
{
  my ($x, $y) = @_;
  
  return (0, 0) if nelem($x) < 2 || nelem($y) < 2 || nelem($x) != nelem($y);
  my $n = nelem($x);
  my $Sxx = sum($x * $x);
  my $Syy = sum($y * $y);
  my $Sxy = sum($x * $y);
  return (0, 0) if $Sxx == 0;
  
  my $a = $Sxy / $Sxx;

  my $SDr = sqrt(($Syy - $a*$Sxy) / ($n - 1));
  my $SEa = $SDr / sqrt($Sxx);

  return $a unless wantarray;
  return $a, $SEa;
}

######################################################

sub wLinFit
#
# Linear weighted fit y = ax + b
#
{
  my ($x, $y, $w) = @_;
  return (0, 0) if nelem($x) < 3 || nelem($y) < 3 || nelem($x) != nelem($y);
  return LinFit($x, $y) unless defined $w;

  my $Sw = sum($w);
  my $Sx = sum($w * $x);
  my $Sy = sum($w * $y);
  my $Sxx = sum($w * $x * $x);
  my $Sxy = sum($w * $x * $y);
  my $delta = $Sw * $Sxx - $Sx**2;
  my $a = ($Sw * $Sxy - $Sx * $Sy) / $delta;
  my $b = ($Sxx * $Sy - $Sx * $Sxy) / $delta;
  
  my $SEa = sqrt($Sw / $delta);
  my $SEb = sqrt($Sxx / $delta);
  
  return ($a, $b, $SEa, $SEb);
}

######################################################

=item corr

  my $r = corr($x, $y);
  my ($r, $p) = corr($x, $y);

Pearson's correlation coefficient (and its significance). $x and $y are piddles of the same size.

=cut
  
sub corr
#
# Pearson's correlation coefficient
#
{
  my ($x, $y, $default) = @_;
  $default = 0 unless defined $default;
  
  my $n = nelem($x);
  return $default if $n < 2;
  my ($mx, $sx) = stats($x);
  my ($my, $sy) = stats($y);
  return $default if $sx == 0 || $sy == 0;
  my $s = sum(($x - $mx) * ($y - $my)) / (nelem($x) - 1);
  my $r =  $s / ($sx * $sy);
  
  my $t = $r / sqrt((1 - $r*$r) / ($n - 2));
  my $p = 1 - pt($t, $n - 2);
  
  return $r unless wantarray;
  return $r, $p;
}

sub corrb
#
# Pearson's correlation coefficient allowing for bad values
#
{
  my ($x, $y, $default) = @_;
  $default = 0 unless defined $default;
  
  my $sel = which(($x->isgood()) & ($y->isgood()));
  $x = $x($sel);
  $y = $y($sel);
  return $default if nelem($x) < 2;
  my ($mx, $sx) = stats($x);
  my ($my, $sy) = stats($y);
  
  $sx=at($sx,0); $sy=at($sy,0);   # temporary solution, there is probably a bug
  #print "mx=$mx  my=$my  sx=$sx  sy=$sy\n";
  
  return $default if $sx == 0 || $sy == 0;
  my $s = sum(($x - $mx) * ($y - $my)) / (nelem($x) - 1);
  #print "### S= $s\n";
  return $s / ($sx * $sy);
}

################################################

=item ucorr

  my $r = ucorr($x, $y);

Unbiased Pearson's correlation coefficient estimator. $x and $y are piddles of the same size. See Olkin & Pratt (1958), Ann. Math. Stat, 29, 201 

=cut

sub ucorr
{
  my ($x, $y, $default) = @_;
  $default = 0 unless defined $default;
    
  my $n = nelem($x);
  return $default if nelem($x) < 2;
  my $r = corr($x, $y, $default);
  
  my @C = (1, 0.25, 0.28125, 0.585938, 1.79443, 7.2675, 36.6401, 221.149, 1554.95,
12482.8, 112658, 1.12914e6, 1.2444e7, 1.49568e8, 1.94705e9, 2.72911e10, 4.09793e11);
  my $F = 1;
  my $m = ($n - 1) / 2;
  my $q = 1 - $r**2;
  for my $k (1 .. @C - 1)
    {$F += $C[$k] * ff($m, $k - 1) * $q**$k}
  return $r * $F 
}

sub ff
{
  my ($m, $i) = @_;
  
  my $f = 1;
  $f *= 1 / ($m + $_) for (0 .. $i);
  return $f
}

################################################

=item Spearman

  my $S = Spearman($x, $y)
  
Spearman's rank correlation coefficient.

=cut

sub Spearman
{
  my ($x, $y, $default) = @_;
  $default = 0 unless defined $default;
  
  return $default if nelem($x) < 2;
  my $r1 = rank($x);
  my $r2 = rank($y);
  return corr($r1, $r2)  
}

sub rank
{
  my $x = shift;
  my $i = qsorti $x;
  my $xs = $x($i);
  my $s = sequence($x);
  
  my $n = nelem($x);
  my $j = 0;
  my @tie = ();  # tie groups
  for my $k (1 .. $n - 1)
  {
    if(at($xs, $k) != at($xs, $k - 1))
    {
      if($k - 1 > $j)   # tied rank
      {
        $s($j:$k-1) .= ($j + $k - 1) / 2;
        push @tie, $k - $j;     # number of ties in this group
      }
      $j = $k;
    }
  }
  if($n - 1 > $j)
    {$s($j:$n-1) .= ($j + $n - 1) / 2}
  $s++;
    
  my $r = zeroes($x);
  $r($i) .= $s;
  return $r unless wantarray;
  return $r, pdl(@tie);
}

sub rank_
# this is fast, but does not take ties into account!
{
  my $x = shift;
  my $i = qsorti $x;
  my $s = sequence($x) + 1;
  my $r = zeroes($x);
  $r($i) .= $s;
  return $r
}

################################################

BEGIN {

my $oldm = -1;         # 'static' variables
my ($sq, $alxm, $g);

=item genPoisson

  my $p = genPoisson($mean);

Generates a number with a Poissonian distribution of a given mean. From Numerical Recipes.

NOTE: a much better generator is available from Math::Random

=cut

sub genPoisson
{
  my $xm = shift;
  
  my $em;
  
  if($xm < 12)
  {
    if($xm != $oldm)
    {
      $oldm = $xm;
      $g = exp(-$xm);
    }
    $em = -1;
    my $t = 1;
    do {
      $em++;
      $t *= rand();
    } while($t > $g);
  }
  
  else  # $xm >= 12
  
  {
    if($xm != $oldm)
    {
      $oldm = $xm;
      $sq = sqrt(2 * $xm);
      $alxm = log($xm);
      $g = $xm * $alxm - Gammaln($xm + 1);
    }
    my $t;
    do {
      my $y;
      do {
        $y = tan($PI * rand());
        $em = $sq * $y + $xm;
      } while($em < 0);
      $em = int($em);
      $t = 0.9 * (1 + $y*$y) * exp($em * $alxm - Gammaln($em + 1) - $g);
    } while(rand() > $t);
  }
  return $em;
}

}

############################################################

sub QKS
{
  my $lambda = shift;

  my $EPS1 = 0.001;
  my $EPS2 = 1e-8;

  my $a2 = -2 * $lambda**2;
  my $fac = 2;
  my $probks = 0;
  my $termbf = 0;
  for my $j (1 .. 100)
  {
    my $term = $fac * exp($a2 * $j**2);
    $probks += $term;
    return $probks if abs($term) <= $EPS1 * $termbf || abs($term) <= $EPS2 * $probks;
    $fac *= -1;
    $termbf = abs($term);
  }
  return 1;
}

############################################################

=item KSone

  my ($D, $Prob) = KSone($data, $f, $m, $s);

Performs Kolmogorov-Smirnov test for one samples and a given probability distribution. Adapted from Numerical Recipes. See chapter 14 for details. Above, I<$data> is a piddle, $f is a function of three parameters, f(x, m, s). This function should return cumulative probability distribution f(x), for the given mean, m, and standard deviation, s. I<$D> is the distance between the samples and I<$Prob> is the probability of distance being greater than the measured value. The null hypothesis is that the samples are drawn from the population defined by f.

=cut

sub KSone
{
  my ($data, $f, $m, $s) = @_;
  
  my $n = nelem($data);
  my $sdata = qsort $data;
  my $D = 0;
  my $fo = 0;
  my ($x, $y1, $y2);
  
  for my $j (1 .. $n)
  {
    my $fn = $j / $n;
    my $ff = &$f(at($sdata, $j-1), $m, $s);
    my $dtn = abs($fn - $ff);
    my $dto = abs($fo - $ff);
    my $dt = ($dtn > $dto) ? $dtn : $dto;
    if($dt > $D)
    {
      $D = $dt;
      $x = at($sdata, $j-1);
      $y1 = $ff;
      $y2 = $fn;
    }
    $fo = $fn;
  }
  
  my $sn = sqrt($n);
  my $Prob = QKS(($sn + 0.12 + 0.11/$sn) * $D);
  return ($D, $Prob, $x, $y1, $y2)
}

############################################################

=item KStwo

  my ($D, $Prob) = KStwo($data1, $data2);

Performs Kolmogorov-Smirnov test for two samples. Adapted from Numerical Recipes. See chapter 14 for details. Above, I<$data1> and I<$data2> are piddles, I<$D> is the distance between the samples and I<$Prob> is the probability of distance being greater than the measured value. The null hypothesis is that both samples are drawn from the same distribution. If KStwo returns a small value of I<$Prob> then one can reject the hypothesis, whith probability I<$Prob> of being wrong. Otherwise, nothing is granted. In particular, a large value of I<$Prob> doesn't validate the hypothesis!

=cut

sub KStwo
# Kolmogorov-Smirnov test
# From numerical recipes
# KS test for two samples
# ($D, $p) = KSTwo($data1, $data2)
# returns distance $D and probability $p
# $data1 and $data2 are piddles
{
  my ($data1, $data2) = @_;
  
  my $n1 = nelem($data1);
  my $n2 = nelem($data2);
  die "\nNeed at least 4 elements in KStwo\n" if $n1 < 4 || $n2 < 4;
  
  my $sdata1 = qsort $data1;
  my $sdata2 = qsort $data2;
  
  my $j1 = 0;
  my $j2 = 0;
  my $fn1 = 0;
  my $fn2 = 0;
  my $D = 0;
  my ($x, $y1, $y2);
  
  while($j1 < $n1 && $j2 < $n2)
  {
    my $d1 = at($sdata1, $j1);
    my $d2 = at($sdata2, $j2);

    $fn1 = ++$j1 / $n1 if $d1 <= $d2;
    $fn2 = ++$j2 / $n2 if $d2 <= $d1;
    
    my $dt = abs($fn2 - $fn1);
    if($dt > $D)
    {
      $D = $dt;
      $x = ($d1 + $d2) / 2;
      $y1 = $fn1;
      $y2 = $fn2;
    }
  }
  my $sn = sqrt($n1 * $n2 / ($n1 + $n2));
  my $Prob = QKS(($sn + 0.12 + 0.11/$sn) * $D);
  
  return ($D, $Prob, $x, $y1, $y2);
}

###############################################

=item KSone2D

  my ($D, $Prob) = KSone2D($x, $y, $f);

Performs two-dimentional Kolmogorov-Smirnov test for one sample and a given probability distribution. Adapted from Numerical Recipes. See chapter 14 for details.

$x, $y - input 2D sample

$f - a subroutine that returns a piddle (fa, fb, fc, fd) with fractions of the cumulative distribution from the four quadrants. The quadrants are defined as

  b | a
  -----
  c | d 

The quadrant function for the uniform distribution over the unit square is provided as C<quaduni1>. So, to perform the 2D K-S test against the uniform distribution you need to call 

  my ($D, $Prob) = KSone2D($x, $y, \&quaduni1);
  
B<Warning:> this test is not very accurate for n < 20. Also, if the returned p-value is > 0.2, it is not accurate.
  
This 2D K-S test is based on
  - Fasano G. & Franceschini A. (1987), MNRAS, 225, 155
  - Peacock J. A. (1983), MNRAS, 202, 615

=cut

sub KSone2D
{
  my ($xx, $yy, $f) = @_;
  
  my $n = nelem($xx);
  die "Unequal number of elements in x and y in KSone2D\n" unless nelem($yy) == $n;
  
  my $D = 0;
  for my $i (0 .. $n - 1)
  {
    my $x = at($xx, $i);
    my $y = at($yy, $i);
    
    my $ff = quadct($x, $y, $xx, $yy);
    my $gg = &$f($x, $y);
    my $m = max(abs($ff - $gg));
    $D = $m if $m > $D;
  }
  my $r = corr($xx, $yy);
  my $sn = sqrt($n);
  my $Prob = QKS($D * $sn / (1 + sqrt(1 - $r**2) * (0.25 - 0.75/$sn)));
  
  return $D, $Prob
}

##############################################

# NOT TESTED!

sub KStwo2D
{
  my ($xx1, $yy1, $xx2, $yy2) = @_;
  
  my $n1 = nelem($xx1);
  die "Unequal number of elements in x1 and y1 in KStwo2D\n" unless nelem($yy1) == $n1;
  my $n2 = nelem($xx2);
  die "Unequal number of elements in x2 and y2 in KStwo2D\n" unless nelem($yy2) == $n2;
  
  my $D1 = 0;
  for my $i (0 .. $n1 - 1)
  {
    my $x = at($xx1, $i);
    my $y = at($yy1, $i);
    
    my $f1 = quadct($x, $y, $xx1, $yy1);
    my $f2 = quadct($x, $y, $xx2, $yy2);
    my $m = max(abs($f1 - $f2));
    $D1 = $m if $m > $D1;
  }
  
  my $D2 = 0;
  for my $i (0 .. $n2 - 1)
  {
    my $x = at($xx2, $i);
    my $y = at($yy2, $i);
    
    my $f1 = quadct($x, $y, $xx1, $yy1);
    my $f2 = quadct($x, $y, $xx2, $yy2);
    my $m = max(abs($f1 - $f2));
    $D2 = $m if $m > $D2;
  }
  
  my $D = ($D1 + $D2) / 2;
  my $r1 = corr($xx1, $yy1);
  my $r2 = corr($xx2, $yy2);
  my $rr = sqrt(1 - ($r1**2+$r2**2)/2);
  my $sn = sqrt($n1 * $n2 / ($n1 + $n2));
  my $Prob = QKS($D * $sn / (1 + $rr * (0.25 - 0.75/$sn)));
  
  return $D, $Prob
}

##############################################

# Quadrants
#
#  b | a
# -------
#  c | d 
#

sub quadct
# Part of 2D Kolmogorov-Smirnof test
# from Numerical Recipes
# returns the fraction of points in each quadrant
{
  my ($x, $y, $xx, $yy) = @_;
  
  my $n = nelem($xx);
  die "Unequal number of elements in xx and yy in quadct\n" unless nelem($yy) == $n;

  my $na = nelem(which(($yy > $y) & ($xx > $x)));
  my $nb = nelem(which(($yy > $y) & ($xx <= $x)));
  my $nc = nelem(which(($yy <= $y) & ($xx <= $x)));
  my $nd = nelem(which(($yy <= $y) & ($xx > $x)));
  
  return pdl($na/$n, $nb/$n, $nc/$n, $nd/$n)
}

sub quaduni1
# quadrant fractions for a unit square
# returns the fraction of area in each quadrant
{
  my ($x, $y) = @_;
  
  my $fa = (1 - $x) * (1 - $y);
  my $fb = $x * (1 - $y);
  my $fc = $x * $y;
  my $fd = (1 - $x) * $y;
  
  return pdl($fa, $fb, $fc, $fd)
}

##############################################

=item tTest

  my $p = tTest($x1, $x2);
  my ($p, $t, $dof) = logTTest($x1, $x2);
  
Standard t-test for two samples with identical variances.

=cut

sub tTest
{
  my ($d1, $d2) = @_;
  
  my $mean1 = mean($d1);
  my $mean2 = mean($d2);
  my $n1 = nelem($d1);
  my $n2 = nelem($d2);
  my $S1 = sum(($d1 - $mean1)**2);
  my $S2 = sum(($d2 - $mean2)**2);
 
  my $sd = sqrt(($S1 + $S2) / ($n1 + $n2 - 2) * (1/$n1 + 1/$n2));
  my $t = abs($mean1 - $mean2) / $sd;
  my $dof = $n1 + $n2 - 2;
  #my $p = Betai(0.5*$dof, 0.5, $dof/($dof + $t**2));
  my $p = 2*(1 - pt($t, $dof));
  return $p unless wantarray;
  return ($p, $t, $dof);
}

#########################################################

=item logtTest

  my $p = logTTest($x1, $x2);
  my ($p, $t, $dof) = tuTest($x1, $x2);
  
Log-ratio t-test appropriate for log-normal distribution.

=cut

sub logtTest
{
  my ($x1, $x2) = @_;
  
  my $n1 = nelem($x1);
  my $n2 = nelem($x2);
  my ($m1, $s1) = stats($x1);
  my ($m2, $s2) = stats($x2);
  
  my $t = abs(log($m2) - log($m1)) / sqrt((($s1/$m1)**2)/$n1 + (($s2/$m2)**2)/$n2);
  my $dof = $n1 + $n2 - 2;
  my $p = 2*(1 - pt($t, $dof));
  return $p unless wantarray;
  return ($p, $t, $dof);
}

# this is for testing only

sub logntTest
{
  my ($x1, $x2) = @_;
  
  my $n1 = nelem($x1);
  my $n2 = nelem($x2);
  my ($m1, $s1) = stats($x1);
  my ($m2, $s2) = stats($x2);
  
  my $t = abs(log($m2) - log($m1)) / sqrt((($s1/$m1)**2)/$n1 + (($s2/$m2)**2)/$n2);
  my $dof = $n1 + $n2 - 2;
  #my $p = 2*(1 - pt($t, $dof));
  my $p = 2*(1 - pnorm($t));
  return $p unless wantarray;
  return ($p, $t, $dof);
}
#########################################################

=item OneSampletTest

  my $p = OneSampletTest($x, $x0);
  my ($p, $t, $dof) = OneSampletTest($x, $x0);
  
One sample t-test. The null hypothesis is that mean($x) = $x0.

=cut

sub OneSampletTest
{
  my ($d, $d0) = @_;
  $d0 ||= 0;
  
  my ($m, $s) = stats($d);
  my $n = nelem($d);
  my $t = abs($m - $d0) / ($s / sqrt($n));
  my $dof = $n - 1;
  my $p = 2*(1 - pt($t, $dof));
  return $p unless wantarray;
  return ($p, $t, $dof);
  
}

#########################################################

=item tuTest

  my $p = tuTest($x1, $x2);
  my ($p, $t, $dof) = tuTest($x1, $x2);
  
Standard t-test for two samples with different variances.

=cut

sub tuTest
{
  my ($d1, $d2) = @_;
  
  my ($mean1, $stdev1) = stats($d1);
  my ($mean2, $stdev2) = stats($d2);
  my $var1 = $stdev1**2;
  my $var2 = $stdev2**2;
  my $n1 = nelem($d1);
  my $n2 = nelem($d2);
  
  my $t = abs($mean1 - $mean2) / sqrt($var1/$n1 + $var2/$n2);
  my $dof = ($var1/$n1 + $var2/$n2)**2 / (($var1/$n1)**2/($n1-1) + ($var2/$n2)**2/($n2-1));
  #my $p = Betai(0.5*$dof, 0.5, $dof/($dof + $t**2));
  my $p = 2*(1 - pt($t, $dof));
  return $p unless wantarray;
  return ($p, $t, $dof);
}

#########################################################

=item ShrinktTest

  my $p = ShrinktTest($x1, $x2);
  my ($p, $m1, $vs1, $m2, $vs2) = tuTest($x1, $x2);
  
Shrinkage t-test using shrinkage variance of Opgen-Rhein & Strimmer (2007). $x1 and $x2 are 2D piddles with rows containing "genes" and columns containing "replicates". Shrinkage variance is calculated for each "condition" and t-test is performed per each "gene".

Returns either a piddle of p-values (uncorrected) or p-values, means and shrunk variances.

=cut

sub ShrinktTest
{
  my ($x1, $x2) = @_;

  my ($nrep1, $ngen1) = dims($x1);
  my ($nrep2, $ngen2) = dims($x2);
  croak "Number of rows must be equal in ShrinktTest ($ngen1, $ngen2)\n" unless $ngen1 == $ngen2;

  my $n1 = float(sumover(isgood($x1)));
  my $n2 = float(sumover(isgood($x2)));
  croak "Need at least two good data points in each row in ShrinktTest\n" if min($n1) < 2 || min($n2) < 2;

  my ($m1, $s1) = statsover($x1);
  my ($m2, $s2) = statsover($x2);

  my $vs1 = ShrinkVar($x1);
  my $vs2 = ShrinkVar($x2);
  my $ss12 = sqrt(0.5 * ($vs1 + $vs2));
  my $ts = abs($m1 - $m2) / sqrt($vs1 / $n1 + $vs2 / $n2);
  
  my @Ps = ();
  for my $i (0 .. $ngen1 - 1)
  {
    my $t = at($ts, $i);
    my $dof = at($n1, $i) + at($n2, $i) - 2;
    my $p = (isgood($ts($i))) ? 2*(1 - pt($t, $dof)) : 0; 
    push @Ps, $p 
  }
  my $Ps = pdl(@Ps);
  return $Ps unless wantarray;
  return $Ps, $m1, $vs1, $m2, $vs2;
}

sub OneSampleShrinktTest
{
  my ($x, $x0) = @_;
  $x0 ||= 0;

  my ($nrep, $ngen) = dims($x);
  my $n = float(sumover(isgood($x)));
  croak "Need at least two good data point in each row in OneSampleShrinkTest\n" if min($n) < 2;
  my ($m, $s) = statsover($x);

  my $vs = ShrinkVar($x);
  my $ts = abs($m - $x0) / sqrt($vs / $n);
  
  my @Ps = ();
  for my $i (0 .. $ngen - 1)
  {
    my $t = at($ts, $i);
    my $dof = at($n, $i) - 1;
    my $p = (isgood($ts($i))) ? 2*(1 - pt($t, $dof)) : 0; 
    push @Ps, $p 
  }
  my $Ps = pdl(@Ps);
  return $Ps unless wantarray;
  return $Ps, $m, $vs;
}

=item ShrinkVar

  my $v = ShrinkVar($x);
  
Shrinkage variance of Opgen-Rhein & Strimmer (2007). Input data $x is a 2D piddle, where rows represent "genes" and columns represent "replicates". The shrinkage target is the median of all gene variances.

The output is a 1D piddle with calculated variances.

ShrinkVar does handle bad values.   

R. Opgen-Rhein and K. Strimmer. 2007. Accurate ranking of differentially expressed genes by a distribution-free shrinkage approach. Statist. Appl. Genet. Mol. Biol. 6: 9

=cut

sub ShrinkVar
{
  my ($x) = @_;
  
  my ($N, $G) = dims $x;

  # actual number of "good" replicates per gene
  # float is needed to avoid integer rounding later
  my $n = float(sumover(isgood($x)));
  #croak "Need at least two good data points in each row\n" if min($n) < 2;
  croak "Need at least one good data point in each row\n" if min($n) < 1;

  my $xm = transpose(sumover($x) / $n);  # \bar{x}_k
  my $w = ($x - $xm)**2;                           # w_ik
  my $wm = sumover($w) / $n;                  # \bar{w}_k
    
  my $v = zeroes($wm);
  my $sel = which($n > 1);
  $v($sel) .=  $n($sel) / ($n($sel) - 1) * $wm($sel);   # v_k
  #my $v = $n / ($n - 1) * $wm;
  
  my $vmed = median($v($sel));
  $wm = transpose($wm);
  my $vt = transpose($v);

  # $var is a 2D piddle with dims 1xN (one column)
  # necessary to do sumover($w(,sel) - $wm(,$sel))
  my $var = zeroes($wm);
  $var(,$sel) .= transpose($n($sel) / ($n($sel) - 1)**3 * sumover(($w(,$sel) - $wm(,$sel))**2));
  #my $var = transpose($n / ($n - 1)**3 * sumover(($w - $wm)**2));
    
  my $lambda = sum($var(,$sel)) / sum(($v($sel) - $vmed)**2);
  $lambda = 1 if $lambda > 1;
  my $vsh = $lambda * $vmed + (1 - $lambda) * $vt;
  #print "lambda = $lambda\n";
  
  #print "x=$x\nn=$n\nxm=$xm\nw=$w\nwm=$wm\nv=$v\nvmed=$vmed\nvar=$var\nlambda=$lambda\nvsh=$vsh\n\n";die;
  
  return $vsh(;-)
}

sub ShrinkVarUnitTest
{
  my $NAN = -666;
  my $x = pdl(
    [1,3,5],
    [3,3,4],
    [2,2,2],
    [1,4,$NAN],
    [1,$NAN,$NAN]
  );

  #expected results:
  my $vr = pdl(3.56205, 0.77128, 0.51757, 3.94261, 0.51757);
  my $eps = 1e-4;
  
  $x->badflag(1);
  $x->badvalue($NAN);
  
  my $v = ShrinkVar($x);
  my $delta = sum(($v - $vr)**2);
  
  return $delta < $eps;
}

#######################################################

=item MannWhitney

  my $p = MannWhitney($x, $y, $noties);

Performs Mann-Whitney (Wilcoxon rank sum) test on samples $x and $y (piddles). Returns p-value from normal distribution. Sample sizes should be greater than 8.

Ties can be taken into account: half a tie is allocated to one group and half to another. An adjusted standard deviation is then calculated.

Sheskin D. J. (2004), Handbook of Parametric and Nonparametric Statistical Procedures: Third Edition, New York: Chapman & Hall, p428-431

=cut

sub MannWhitney
{
  my ($x, $y) = @_;
  
  my $nx = nelem($x);
  my $ny = nelem($y);
  my $n = $nx + $ny;
    
  my $d = $x->append($y);     # pooled sample
  my $lab = zeroes($d);       # labels
  $lab(0:$nx-1) .= 1;         # x -> 1, y -> 0
  
  my ($r, $tie) = rank($d);   # ranks and tied group counts
  my $rx = sum($r->where($lab == 1));
  my $ry = sum($r->where($lab == 0));
  
  my $nxy = $nx * $ny;
  my $ux = $nxy + $nx*($nx + 1)/2 - $rx;
  my $uy = $nxy + $ny*($ny + 1)/2 - $ry;
  my $u = ($ux > $uy) ? $uy : $ux;

  my $M = $nxy / 2;
  my $V = $nxy * ($n + 1) / 12;
  $V = $V * (1 - 1/($n**3 - $n) * sum($tie**3 - $tie));
  
  my $Z = abs($u - $M) / sqrt($V);
  my $p = 2*(1 - pnorm($Z));

  return $p unless wantarray;
  return $p, $u;
}

#######################################################

=item PermutationTest

  my $p = PermutationTest($x, $y, $nsim, $stat, $alpha);

Performs permutation test for two samples. The null hypothesis is that they come from the same distribution. The test is based on $nsim random draws of two samples the same size as $x and $y from pooled data ($x, $y). The fraction of instances when the simulated test statistic is greater than the observed test statistic, is the returned p-value.

There are 4 possible test statistics, $stat:

  C<mean>: D = abs(mean($x) - mean($y))
  C<mean_ratio>: D = abs(ln(mean($x) / mean($y))), mean($y) != 0
  C<median>: D = abs(median($x) - median($y))
  C<tmean>: D = abs(tmean($x, $alpha) - tmean($y, $alpha))
  
where tmean is the trimmed mean with fraction $alpha of data cut off (see TrimmedMean).

=cut

sub PermutationTest
{
  my ($x, $y, $nsim, $stat, $alpha, $verbouse) = @_;
  $nsim ||= 10000;
  $stat ||= 'mean';
  $alpha ||= 0.5;
  
  my $nx = nelem($x);
  my $ny = nelem($y);
  my $n = $nx + $ny;
  my $Dobs = PermutationTestStatistic($x, $y, $stat, $alpha);
  
  my $d = $x->append($y);
  my $cnt = 0;
  for my $i (1 .. $nsim)
  {
    Percent1($i/$nsim) if $verbouse && $i % 1000 == 0;
    my $r = random($n);
    my $idx = qsorti $r;
    my $di = $d($idx);
    my $xi = $di(0:$nx-1);
    my $yi = $di($nx:$n-1);
    my $D = PermutationTestStatistic($xi, $yi, $stat, $alpha);
    $cnt++ if $D > $Dobs;
  }
  print "\n" if $verbouse;
  return $cnt / $nsim
}

sub PermutationTestStatistic
{
  my ($x, $y, $stat, $alpha) = @_;
  
  my $S;
  if($stat eq 'mean')
    {$S = abs(mean($x) - mean($y))}
  elsif($stat eq 'mean_ratio')
  {
    my $mx = mean($x);
    my $my = mean($y);
    # this is for RNASeq data only:
    $mx = 1e-16 if $mx == 0;
    $my = 1e-16 if $my == 0;
    #die "Cannot divide by zero in PermutationTest\n" if $my == 0;
    $S = abs(log($mx / $my))
  }
  elsif($stat eq 'median')
    {$S = abs(median($x) - median($y))}
  elsif($stat eq 'tmean')
    {$S = abs(TrimmedMean($x, $alpha) - TrimmedMean($y, $alpha))}
  else
    {die "Unrecognized test statistic $stat in PermutationTest\n"}
  return $S
}

#######################################################

=item TrimmedMean

  my $m = TrimmedMean($x, $alpha);

Trimmed mean of piddle $x. Mean calculated after discarded fraction $alpha of the sample at high and low end. Data are sorted and $alpha/2 of the points are discarded on each end. If $alpha/2 is not an integer, result is interpolated betwteen the two means of integer truncations.

=cut

sub TrimmedMean
{
  my ($d, $alpha) = @_;
  
  my $low = pct($d, $alpha / 2);
  my $up = pct($d, 1 - $alpha / 2);
  my $tmean = average($d->where(($d >= $low) & ($d <= $up)));
  
  return $tmean; 
}


sub TrimmedMean_
{
  my ($x, $alpha) = @_;
  
  my $s = qsort $x;
  my $n = nelem($s);
  my $ntrim = $n * $alpha / 2;
  my $n1 = floor($ntrim);
  
  my $tmean;
  if($ntrim == $n1)
  {
    $tmean = mean($s($ntrim:$n-$ntrim-1)); 
  }
  else  # non-integer trim, need to interpolate
  {
    my $n2 = $n1 + 1;
    #my $tm1 = mean($s($n1:$n-$n1-1));
    my $tm2 = mean($s($n2:$n-$n2-1));
    # it is faster to add two points than calculate another mean:
    my $tm1 = ($tm2*($n-2*$n2) + at($s,$n1) + at($s,$n-$n1-1)) / ($n-2*$n1);
    $tmean = $tm1 + ($ntrim - $n1) * ($tm2 - $tm1);  # interpolation
  }
  return $tmean
}

sub mean
{
  my $x = shift;
  sum($x) / nelem($x)
}

#######################################################

=item BootstrapTest

  my $p = PermutationTest($x, $y, $B);

Performs bootstrap test for sample means. The test is based on $B random draws with replacement of two samples the same size as $x and $y from pooled data ($x, $y). t statstic is calculated for each bootstrap (t*) and p-value calculated as the fraction of bootstrap t* > t_obs. Based on Efron B., Tibshirani R. J. (1993), An introdoction to the bootstrap, Chapman & Hall. 

=cut

sub BootstrapTest
{
  my ($x1, $x2, $B, $samevar) = @_;
  $B ||= 10000;
  
  # sample statistics
  my ($tobs, $n1, $n2, $M1, $M2) = tb($x1, $x2, $samevar);
  my $n = $n1 + $n2;
  
  # pooled normalized sample
  my $M = ($n1 * $M1 + $n2 * $M2) / $n;   # common mean
  my $T1 = $x1 - $M1 + $M;
  my $T2 = $x2 - $M2 + $M;
  my $T = append($T1, $T2);   # pooled normalized sample
  
  # bootstrap
  my $crit = 0;
  for my $b (1 .. $B)
  {
    # resampling with replacement
    my $Z1 = $T(floor(random($n1) * $n));
    my $Z2 = $T(floor(random($n2) * $n));
    my $tb = tb($Z1, $Z2, $samevar);
   
    $crit++ if $tb >= $tobs;
  }
  my $p = $crit / $B;
  
  return $p unless wantarray; 
  return $p, $M1, $M2;
}

sub tb
{
  my ($x1, $x2, $samevar) = @_;
  
  my $n1 = nelem($x1);
  my $n2 = nelem($x2);
  my $M1 = sum($x1) / $n1;
  my $M2 = sum($x2) / $n2;
  my $S1 = sum(($x1 - $M1)**2);
  my $S2 = sum(($x2 - $M2)**2);

  my $tb;  
  if($samevar)
  {
    my $V12 = ($S1 + $S2) / ($n1 + $n2 - 2);
    $tb = abs($M2 - $M1) / sqrt($V12  * (1/$n1 + 1/$n2));
  }
  else
  {
    my $V1 = $S1 / ($n1 - 1);
    my $V2 = $S2 / ($n2 - 1);
    $tb = abs($M2 - $M1) / sqrt($V1/$n1 + $V2/$n2);
  }
  return $tb unless wantarray;
  return $tb, $n1, $n2, $M1, $M2, $S1, $S2; 
}

###############################################

sub logfact
{
   #return Gammaln(shift(@_) + 1.0);
   CompBio::SpecFun::Gamma::gammaln(shift(@_) + 1.0)
}

#########################################################

sub hypergeom
{
   # There are m "bad" and n "good" balls in an urn.
   # Pick N of them. The probability of i or more successful selections:
   # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
   my ($n, $m, $N, $i) = @_;

   my $loghyp1 = logfact($m)+logfact($n)+logfact($N)+logfact($m+$n-$N);
   my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact($N-$i)+logfact($m+$n);
   return exp($loghyp1 - $loghyp2);
}

#########################################################

sub binomial
{
  my ($n, $k) = @_;
  
  my $logbin = logfact($n) - logfact($k) - logfact($n - $k);
  return int(exp($logbin) + 0.5);
}

#########################################################

=item LogrankTest

  my $p = LogrankTest($n1, $n2)
  
Logrank test for B<uncensored> survival data. $n1 and $n2 are ordered surviving numbers in two groups. Both data are on the same time grid.

=cut

sub LogrankTest
{
  my ($n1, $n2) = @_;
  
  my $N = nelem($n1);
  die "Both data sets need to have the same size in LogrankTest\n" unless nelem($n2) == $N;
  
  my $m1 = $n1(0:$N-2) - $n1(1:$N-1);
  my $m2 = $n2(0:$N-2) - $n2(1:$N-1);
  $m1 = append($m1, 0);
  $m2 = append($m2, 0);
  my $sel = which(($m1 > 0) | ($m2 > 0));
  
  #$m1 = pdl(0,0,0,0,0,3,1,0,1,0,0,1,0,1,0,1,1);
  
  $n1 = $n1($sel);
  $n2 = $n2($sel);
  my $n = $n1 + $n2;
  $m1 = $m1($sel);
  $m2 = $m2($sel);
  my $m = $m1 + $m2;
  
  my $e1 = $m * $n1 / $n;    # expected no of events in 1
  my $me1 = $m1 - $e1;       # O - E
  my $v1 = $e1 * (1 - $n1 / $n) * ($n - $m) / ($n - 1);   # variance 
  my $OE = sum($me1);
  my $V = sum($v1);
  my $X = $OE**2 / $V;  # chi-square distributed variable
  
  #my $E = sum($e1);
  #my $Y = $OE**2 / $E; 
  
  my $p = 1 - pchisq($X, 1);
  
  #wcols $m1, $m2, $n1, $n2, $me1;
  #print "E = $E\nOE = $OE\nOEO = $Y\nX2 = $X\n";
  #print "p = $p\n";
  return $p
}


#########################################


=item EnrichmentTest

  my $p = EnrichmentTest($nsel, $nuni, $Nsel, $Nuni, $lf)
  
Performs hypergeometric test for enrichment ,e.g., for genes with GO terms.

  $nsel - number of genes in the selected set annotated with a GO term
  $nuni - number of genes in the universe annotated with a GO term
  $Nsel - size of the selected set
  $Nuni - size of the gene universe
  $lf - precalculated log factorial, using InitLogfact(); if not defined, they will be calculated directly

=cut

sub EnrichmentTest
{
  my ($nsel, $nuni, $Nsel, $Nuni, $lf) = @_;
  my $expect = $nuni * $Nsel / $Nuni;
  my $p = 1;
  
  if($nsel > 0)
  {
    my $sum = 0;
    my $min = $nsel;
    my $max = ($nuni < $Nsel) ? $nuni : $Nsel;
    for my $k ($min .. $max)
      {$sum += fast_hypergeom($nuni, $Nuni - $nuni, $Nsel, $k, $lf)};
    $p = $sum
  }
  
  return $p unless wantarray;
  return $p, $expect
}

sub UnitTestEnrichmentTest
{
  my $lf = InitLogfact();
  my $pass = 1;
  # tested against GraphPad online Fisher's test
  my ($p, $exp) = EnrichmentTest(7, 300, 10, 1000, $lf);
  $pass = 0 if abs($p - 0.0102) > 0.0001;
  ($p, $exp) = EnrichmentTest(35, 600, 50, 1000, $lf);
  $pass = 0 if abs($p - 0.0898) > 0.0001;
  ($p, $exp) = EnrichmentTest(45, 400, 100, 1000, $lf);
  $pass = 0 if abs($p - 0.1664) > 0.0001;
  
  if($pass) {print "Passed\n"} else {print "Failed\n"}
  return $pass
}

sub InitLogfact
{
  my ($max) = @_;
  $max ||= 10000;
  
  my @s = ();
  for my $i (0 .. $max)
    {push @s, logfact($i)}
  return \@s;
}

sub fast_hypergeom
{
   # There are m "bad" and n "good" balls in an urn.
   # Pick N of them. The probability of i or more successful selections:
   # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
   my ($n, $m, $N, $i, $logfact) = @_;
   
   return hypergeom($n, $m, $N, $i) unless defined $logfact;

   my $loghyp1 = $logfact->[$m]+$logfact->[$n]+$logfact->[$N]+$logfact->[$m+$n-$N];
   my $loghyp2 = $logfact->[$i]+$logfact->[$n-$i]+$logfact->[$m+$i-$N]+$logfact->[$N-$i]+$logfact->[$m+$n];
   return exp($loghyp1 - $loghyp2);
}


#########################################

=item MedianCI

  my ($median, $ml, $mu) = MedianCI($x, $gamma)

Calculates sample median and its confidence limits, using the method of Hettmansperger T.P. and Sheather S.J. 1986, Stat. Prob. Let. 4, 75. I<$x> is the sample (piddle) and I<$gamma> is the confidence probability (e.g. 0.95).

=cut

sub MedianCI
{
  my ($x, $gamma) = @_;
  
  $gamma ||= 0.95;
  
  my $n = nelem($x);
  my ($mean, $stdev, $median) = stats($x);
  my $xs = qsort($x);
  
  my $d = ceil($n / 2) - 1;
  my $g1 = gammad($d + 1, $n);
  my $g = gammad($d, $n);
  croak "Gamma too small, cannot find confidence limits\n" if $gamma < $g1;
  while($gamma > $g)
  {
    $d--;
    $g1 = $g;                         # gamma_{d+1}
    $g = gammad($d, $n);              # gamma_{d}
  }

  my ($ml, $mu);
  if($d == 0)
  {
    carp "Warning: median confidence limits outside the sample (sample too small)\n";
    $ml = at($xs, 0);
    $mu = at($xs, $n - 1);
  }
  else
  {
    my $I = ($g - $gamma) / ($g - $g1);
    my $lambda = (($n - $d) * $I) / ($d + ($n - 2*$d)*$I);

    my $xd = at($xs, $d - 1);
    my $xd1 = at($xs, $d + 1 - 1);
    my $xnd = at($xs, $n - $d - 1);
    my $xnd1 = at($xs, $n - $d + 1 - 1);
    $ml = (1 - $lambda) * $xd + $lambda * $xd1;
    $mu = (1 - $lambda) * $xnd1 + $lambda *$xnd;
  }
  
  return $median, $ml, $mu;
}

sub gammad
#
# gammad = P(d <= W < n - d), where W is binomial (n, p = 0.5)
#
{
  my ($d, $n) = @_;

  # pbinom returns P(X <= d)
  ($d == 0) ? 1 : 1 - 2*pbinom($d - 1, $n, 0.5)
}

#############################################

=item PropCI

  my ($p, $w) = PropCi($n, $s, $alpha)

Calculates confidence interval of a proportion.

  $n - total number
  $s - number of "successess"
  $alpha - confidence limit (default 0.05)
  
Returns

  $p - modified proportion
  $w - width of the CI: [$p - $w, $p + $w]

=cut


sub PropCI
{
  my ($n, $s, $alpha) = @_;
  $alpha ||= 0.05;
  
  my $Z = qnorm(1 - $alpha / 2);
  my $p = ($s + $Z) / ($n + $Z*$Z);
  $p = 1 if $p > 1;
  $p = 0 if $p < 0;
  my $w = $Z * sqrt($p * (1 - $p) / ($n + $Z*$Z));

  return $p, $w  
}

###########################################

=item CohensKappa

  my ($k, $SE) = CohensKappa($d1, $d2)

Calculates Cohen's kappa correlation between $d1 and $d2.

  $d1, $d2 - 1D piddles with equal number of elements, containing 0s and 1s.
  
Returns

  $k - Cohen's kappa
  $SE - standard error found assuming normal distribution of kappa, which should be fine if number of agreements > 5.

=cut

sub CohensKappa
{
  my ($d1, $d2) = @_;
  
  my $Nyy = nelem(which($d1 & $d2));
  my $Nyn = nelem(which($d1 & !$d2));
  my $Nny = nelem(which(!$d1 & $d2));
  my $Nnn = nelem(which(!$d1 & !$d2));
  my $N = $Nyy + $Nyn + $Nny + $Nnn;
  
  print "\nYY = $Nyy\nYN = $Nyn\nNY = $Nny\nNN = $Nnn\n";
  
  # probability of agreement
  my $Pa = ($Nyy + $Nnn) / $N;
  
  my $f1y = ($Nyy + $Nyn) / $N;   # fraction of cases where R1 said yes
  my $f2y = ($Nyy + $Nny) / $N;   # fraction of cases where R2 said yes
  my $f1n = ($Nny + $Nnn) / $N;   # fraction of cases where R1 said no
  my $f2n = ($Nyn + $Nnn) / $N;   # fraction of cases where R2 said no
  my $Py = $f1y * $f2y;   # prob. that both say yes randomly
  my $Pn = $f1n * $f2n;   # prob. that both say no randomly
  
  # probability of random agreement
  my $Pe = $Py + $Pn;
  
  #kappa
  my $k = ($Pa - $Pe) / (1 - $Pe);
  my $SE = sqrt($Pa * (1 - $Pa) / ($N * (1 - $Pe)**2));
  
  print "k = $k, SE = $SE\n";
  
  return $k, $SE
  
  #printf "Pa = %5.3f\nPe = %5.3f\n\nKappa = %5.3f\nSE = %5.3f\n", $Pa, $Pe, $k, $SE; 
}

=item Disagreemet

  my $D = Disagreement($x, $y);
  my $D = Disagreement($x, $y, log=>1);
  
=cut
 
sub Disagreement
{
  my ($x, $y, %opt) = @_;
  
  my $sel = which(($x > 0) & ($y > 0));
  my $distance;
  
  if($opt{log})
  {  
    my $rat = log($y($sel) / $x($sel)) / log(2);
    $rat -= avg($rat) if $opt{recentre} || $opt{recenter};
    $distance =  sqrt(avg($rat**2));
  }
  elsif($opt{ds})
  {
     my $rat = $y($sel) / $x($sel);
     $rat->where($rat < 1) .= 1/$rat->where($rat < 1);
     $distance =  avg($rat - 1);
  }
  else
  {  
     my $rat = $y($sel) / $x($sel);
     $rat->where($rat < 1) .= 1/$rat->where($rat < 1);
     $distance =  sqrt(avg(($rat - 1)**2));
  }
  return $distance
}


1;
