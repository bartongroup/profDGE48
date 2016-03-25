package Distribution;
require Exporter;
require File::Temp;
require PDL;
require PDL::NiceSlice;
require PDL::Options;
require Math::CDF;
require File::Temp;
require Carp;
require Tools;
require Stats;
require Beta;
require Gamma;
require Statistics::R;

=head1 NAME

CompBio::Distribution - tools for probability distirbutions

=head1 SYNOPSIS

  use PDL;
  use CompBio::Distribution;
  
  my $x = grandom(100);
  my ($loglik, %pars) = FitDistr($x, 'normal');
  my ($f, $F) = GetRFunc($x, 'normal', %pars);
  my ($M, $S) = GetRMS('negative binomial', %pars);
  
=head1 DESCRIPTION
  

=head1 FUNCTIONS

=over 4

=cut

use strict;
use warnings;

use File::Temp qw(tempfile mktemp);
use PDL;
use PDL::NiceSlice;
use PDL::Options;
use Math::CDF qw(:all);
use File::Temp qw(tempfile mktemp);
use Carp;
use Tools;
use Stats;
use Beta qw(beta betainc);
use Gamma qw(gammaln gamma);
use Statistics::R;

our $VERSION = 1.0;

my $PI = 4 * atan(1);

our @ISA = ("Exporter");
our @EXPORT = qw(
  FitDistr
  GetRFunc
  GetRMS
  gauss
  fgauss
  lnorm
  flnorm
  fpoisson
  nbinom
  nbinomms
  fnbinom
  fnbinomms
  fgamma
  fgammams
  logfact
  lnL
  CreateDiscreteGrid
  CreateNBinomialTables
  GenerateNBinomial
  BenjaminiHochberg
  BenjaminiHochbergLimit
  HolmBonferroni
  HolmBonferroniLimit
  PoissonTest
  NBBootstrapTest
  NBKSTest
  NBMiTest
  MeintanisStatistic
  AndersonDarlingNBStatistic
  HinzGurland
);

################################################################

=item FitDistr

  my ($loglik, %pars) = FitDistr($x, $func, %opt);

Fit probability distribution $func to data $x, using maximum likelihood method. This script calles R function C<fitdistr>. $x is a one-dimentional piddle with data. $func is a string describing probability distribution, e.g. C<negative binomial>, C<lognormal>, C<poisson>, C<gamma>, etc.

It returns log-likelihood and best-fitting parameters in a hash:

  $pars{$parname} = [$value, $error]

=cut

sub FitDistr
{
  my ($x, $func, %opt) = @_;

  my $cx = 'c(' . join(',', list($x)) . ')';

  my $cmdfit = <<EOF;
library("MASS")
x = $cx
out = fitdistr(x, "$func")
EOF

  my $cmdout = <<EOF;
print(out)
logLik(out)
warnings()
EOF

  my $R = Statistics::R->new();
  $R->run($cmdfit);
  my $out = $R->run($cmdout);
  
  #print "\n====\n$out\n======\n";
  
  my @out = split(/\n/, $out);
  return 0, (error=>$out[0]) if $out[0] =~ /Error/;

  my @name = getline(shift @out);
  my @val = getline(shift @out);
  my @err = getline(shift @out);
  my %pars = ();
  for my $i (0 .. @name - 1)
    {$pars{$name[$i]} = [$val[$i], $err[$i]]}

  my $line = shift @out;
  my @s = split " ", $line;
  my $loglik = $s[2];
  
  #print "L=$loglik\n";
  return $loglik, %pars;
  
  sub getline
  {
    my ($l) = @_;
    $l =~ s/^\s+//;
    $l =~ s/[\(\)]//g;
    my @s = split " ", $l
  } 
}

###############################################

=item GetRFunc

  my ($f, $F) = GetRFunc($x, $func, %pars);

Returns probability density, $f, and cumulative density, $F, functions for a probability distribution given by string $func and parameters %pars. $f anf $F are calculated for data points $x. See C<FitDistr> for more details on $func and %pars.

=cut


sub GetRFunc
{
  my ($x, $func, %pars) = @_;
  
  if($func eq 'normal')
  {
    my $m = $pars{mean}[0];
    my $s = $pars{sd}[0];
    my $f = fgauss($x, $m, $s);
    my $F = zeroes($x);
    for my $i (0 .. nelem($x) - 1)
    {
      $F($i) .= gauss(at($x, $i), $m, $s);
    }
    return $f, $F;
  }
  elsif($func eq 'lognormal')
  {
    my $m = $pars{meanlog}[0];
    my $s = $pars{sdlog}[0];
    
    my $s2 = $s*$s;
    my $M = exp($m + $s2/2);
    my $S = sqrt((exp($s2) - 1) * exp(2*$m + $s2));
    #print "LOG M=$M S=$S\n";
    
    my $f = zeroes($x);
    my $F = zeroes($x);
    for my $i (0 .. nelem($x) - 1)
    {
      my $xx = at($x, $i);
      $f($i) .= ($xx>0) ? flnorm($xx, $M, $S) : 0;
      $F($i) .= ($xx>0) ? lnorm($xx, $M, $S) : 0;
    }
    return $f, $F;
  }
  elsif($func eq 'negative binomial')
  {
    my $mu = $pars{mu}[0];
    my $r = $pars{size}[0];
    my $p = $mu / ($r + $mu);
    my $f = zeroes($x);
    my $F = zeroes($x);
    for my $i (0 .. nelem($x) - 1)
    {
      $f($i) .= fnbinom(at($x, $i), $r, $p);
      $F($i) .= nbinom(at($x, $i), $r, $p);
      # $f($i) .= enbinom(at($x, $i), $r, $p);
      # $F($i) .= senbinom(at($x, $i), $r, $p);
    }
    return $f, $F;
  }
  elsif($func eq 'poisson')
  {
    my $mu = $pars{lambda}[0];
    my $f = zeroes($x);
    my $F = zeroes($x);
    for my $i (0 .. nelem($x) - 1)
    {
      $f($i) .= fpoisson(at($x, $i), $mu);
      $F($i) .= ppois(at($x, $i), $mu);
    }
    return $f, $F;
  }
  elsif($func eq 'gamma')
  {
    my $alpha = $pars{shape}[0];
    my $beta = $pars{rate}[0];
    my $f = zeroes($x);
    my $F = zeroes($x);
    for my $i (0 .. nelem($x) - 1)
    {
      $f($i) .= fgamma(at($x, $i), $alpha, $beta);
    }
    return $f, $F;
  }
  else
    {die}
}

#####################################################

=item GetRMS

  my ($M, $S) = GetRMS($func, %pars);

Returns mean, $M, and standard deviation, $S, for a probability distribution given by string $func and parameters %pars. This is a simple conversion of original parameters %pars, into more ($M, $S). $func can be C<normal>, C<lognormal>, C<negative binomial>, C<poisson> or C<gamma>.

=cut

sub GetRMS
{
  my ($func, %pars) = @_;
  
  return (0, 0) if scalar (keys %pars) == 0;
  my ($M, $S);
  if($func eq 'normal')
  {
    $M = $pars{mean}[0];
    $S = $pars{sd}[0];
  }
  elsif($func eq 'lognormal')
  {
    my $m = $pars{meanlog}[0];
    my $s = $pars{sdlog}[0];
    
    my $s2 = $s*$s;
    $M = exp($m + $s2/2);
    $S = sqrt((exp($s2) - 1) * exp(2*$m + $s2));
  }
  elsif($func eq 'negative binomial')
  {
    $M = $pars{mu}[0];
    my $r = $pars{size}[0];
    my $p = $M / ($r + $M);
    $S = sqrt($p * $r) / (1 - $p);
  }
  elsif($func eq 'poisson')
  {
    $M = $pars{lambda}[0];
    $S = sqrt($M);
  }
  elsif($func eq 'gamma')
  {
    my $a = $pars{shape}[0];
    my $s = 1 / $pars{rate}[0];
    $M = $a * $s;
    $S = sqrt($a) * $s;
  }
  
  return $M, $S
}

##################################################

=item lnL

  my ($lnL, $sic) = lnL($x, $func, %pars);

Returns log-likelihood and Schwarz information criterion (SIC) for a given data and distribution.

=cut

sub lnL
{
  my ($x, $func, %pars) = @_;
  
  return (0, 0) if scalar (keys %pars) == 0;
  my $n = nelem($x);
  my $Npars = scalar keys %pars;
  
  my ($f, $F) = GetRFunc($x, $func, %pars);
  my $lnL = sum(log($f));
  my $sic = $Npars * log($n) - 2 * $lnL;
  
  return $lnL, $sic;
}

##################################################

=item gauss

  my $P = gauss($x, $m, $s);

Gaussian CDF.

=cut

sub gauss
{
  my ($x, $m, $s) = @_;
  
  my $z = ($x - $m) / $s;
  pnorm($z);
}

##################################################

=item fgauss

  my $P = fgauss($x, $m, $s);

Gaussian PDF.

=cut

sub fgauss
#
# Gaussian PDF
#
{
  my ($x, $m, $s) = @_;
  
  my $s2 = $s*$s;
  (1 / sqrt(2 * $PI * $s2)) * exp(-($x - $m)**2 / (2 * $s2))
}

##################################################

=item lnorm

  my $P = lnorm($x, $m, $s);

Log-normal CDF.

=cut

sub lnorm
{
  my ($x, $M, $S) = @_;
  
  my $s2 = log($S*$S * exp(-2 * log($M)) + 1);
  my $m = log($M) - $s2 / 2;
  my $s = sqrt($s2);

  pnorm((log($x) - $m) / $s);  
}

##################################################

=item flnorm

  my $P = flnorm($x, $m, $s);

Log-normal PDF.

=cut

sub flnorm
#
# Log-normal PDF
#
{
  my ($x, $M, $S) = @_;
  
  my $s2 = log($S*$S * exp(-2 * log($M)) + 1);
  my $m = log($M) - $s2 / 2;
  
  (1 / ($x * sqrt(2 * $PI * $s2))) * exp(-(log($x) - $m)**2 / (2 * $s2))
}

##################################################

=item fpoisson

  my $P = nbinom($x, $m);

Continuous Poisson PDF.

=cut

sub fpoisson
{
  my ($x, $M) = @_;
  
  if($x < 100)
    {$M**$x * exp(-$M) / gamma($x + 1)}
  else
    {fgauss($x, $M, sqrt($M))}
}

##################################################

=item nbinom

  my $P = nbinom($k, $r, $p);

Approximate exact formula for negative binomial CDF. Uses beta function.

=cut

sub nbinom
#
# Negative binomial CDF
#
{
  my ($k, $r, $p) = @_;
  
  1 - betainc($p, $k + 1, $r);  # negative binomial using Math::SpecialFun::Beta 
}

##################################################


=item nbinomms

  my $P = nbinom($k, $m, $s);

Negative binomial CDF, parameterized by mean and standard deviation.

=cut

sub nbinomms
{
  my ($k, $m, $s) = @_;

  return 0 if $s*$s <= $m;
  #die "S^2 <= M in msbinom\n" if $s*$s <= $m;
  return 0 if $k < -0.99;
  my $q = $s*$s - $m;
  my $r = $m*$m / $q;
  my $p = $q / ($s*$s);
  #print "p=$p r=$r\n";
  #mynbinom($k, $r, $p);
  1 - betainc($p, $k + 1, $r);  # negative binomial using Math::SpecialFun::Beta 
}

##################################################

=item fnbinom

  my $P = fnbinom($k, $r, $p);

Approximate exact formula for negative binomial PDF. Uses gamma function.

=cut

sub fnbinom_
{
  my ($k, $r, $p) = @_;
  my $P = mybinomial($k + $r - 1, $k) * $p**$k * (1 - $p)**$r;
  #print "p=$p   k=$k  r=$r  p^k=", $p**$k, "q^r=", (1-$p)**$r, "\n";
  
  return $P;
}

# more stable

sub fnbinom
{
  my ($k, $r, $p) = @_;
  my $P = ($k > 0) ? nbinom($k, $r, $p) - nbinom($k - 1, $r, $p) :
                     nbinom($k, $r, $p);
    
  return $P;
}

#################################################


=item fnbinomms

  my $P = fnbinomms($k, $m, $s);

Negative binomial PDF, parameterized by mean and standard deviation.

=cut

sub fnbinomms
{
  my ($x, $M, $S) = @_;
  
  my $q = $S*$S - $M;
  my $r = $M*$M / $q;
  my $p = $q / ($S*$S);
  #print "$M $S $x $r $p\n";
  fnbinom($x, $r, $p);
}

#########################################################

sub mybinomial
{
  my ($n, $k) = @_;
  
  my $logbin = logfact($n) - logfact($k) - logfact($n - $k);
  #print $logbin, " ", exp($logbin), "\n";
  #return int(exp($logbin) + 0.5);
  return exp($logbin)
}

#########################################################

sub logfact
{
  return gammaln(shift(@_) + 1.0);
}

##################################################

=item fgamma

  my $P = fgamma($x, $alpha, $beta);

Gamma distribution PDF. 

=cut

sub fgamma
{
  my ($x, $alpha, $beta) = @_;
  
  ($beta**$alpha / gamma($alpha)) * $x**($alpha - 1) * exp(-$beta * $x);
}

##################################################


=item fgammams

  my $P = fgammams($x, $m, $s);

Gamma distribution PDF, parameterized by mean and standard deviation. 

=cut

sub fgammams
{
  my ($x, $M, $S) = @_;
  
  my $beta = $M / ($S*$S);
  my $alpha = $M * $beta;
  ($beta**$alpha / gamma($alpha)) * $x**($alpha - 1) * exp(-$beta * $x);
}

##################################################

sub CreateDiscreteGrid
{
  my ($f, $M, $S, $cumul, $maxn, $plimit) = @_;
  
  $maxn ||= 1000;
  $plimit ||= 1e-4;
  
  $plimit *= &$f($M, $M, $S) unless $cumul;
  
  my $min = $M - $S;
  while($min >=0 &&  &$f($min, $M, $S) > $plimit)
    {$min -= 0.5 * $S}
  $min = 0 if $min < 0;

  my $max = $M + $S;
  while((($cumul) ? 1 - &$f($max, $M, $S) : &$f($max, $M, $S)) > $plimit)
    {$max += 0.5 * $S}

  $min = floor($min);
  $max = ceil($max);
  my $n = $max - $min + 1;
  #print "$min - $max\n";
    
  my $gap = 1;
  if($n > $maxn)
    {$gap = int($n / $maxn)}
  
  my $x = sequence(int($n / $gap) + 1) * $gap + $min;
  
  return $x
}

#########################################################

sub CreateNBinomialTables
{
  my ($M, $S) = @_;
  
  carp "Requires S > sqrt(M)\n" if $S <= sqrt($M);
  
  my $x = CreateDiscreteGrid(\&nbinomms, $M, $S, 1);
  my $N = nelem($x);
  my $y = zeroes($x);
  for my $i (0 .. $N - 1)
    {$y($i) .= nbinomms(at($x, $i), $M, $S)}
  #wcols $x, $y;
    
  return $x, $y
}

##################################################

sub GenerateNBinomial
#
# This is my subroutine, using tabulated cumulative distribution
#
{
  my ($x, $y, $n) = @_;
  
  my $r = random($n);
  
  my $b = PDL->null;
  my $err = PDL->null;
  interpolate($r, $y, $x+0.5, $b, $err);
  $b = floor($b + 0.5);
  $b->where($b<0) .= 0;   # side effect of interpolation
  #my $i = vsearch($r, $y);
  #my $b = $x($i);
  return $b;
}

##################################################

sub CreateNBinomialTables_
{
  my ($M, $S, $width) = @_;
  
  $width ||= 10;
  
  #my $N = 1000;
  #my $x = ($M + 6*$S)*sequence($N)/$N;
  #print "$M, $S\n";
  my $x = sequence(int($M + $width*$S));
  #print "x=$x\n";
  my $N = nelem($x);
  my $y = zeroes($x);
  for my $i (0 .. $N - 1)
    {$y($i) .= nbinomms(at($x, $i), $M, $S)}
  return $x, $y
}

##################################################

sub GenerateNBinomial_
{
  my ($x, $y, $n) = @_;
  
  my $r = random($n);
  #my $b = ceil(interpol($r, $y, $x) + 0.5);
  my $i = vsearch($r, $y);
  my $b = $x($i)->copy();
}

##################################################

=item BenjaminiHochberg

  my $pc = BenjaminiHochberg($p);

Corrects a list of p-values (piddle $p) for multiple testing using step-up Benjamini-Hochberg correction. Returns a piddle with corrected p-values.

=cut

sub BenjaminiHochberg
{
  my $p_ = shift;
  my $p = $p_->copy();
  
  my $n = nelem($p);
  my $idx = qsorti $p;
  my $c = $n / (sequence($n) + 1);
  #print "$c\n";
  $p($idx) *= $c;
  for my $j (0 .. $n - 2)
  {
    my $i = at($idx, $n - 2 - $j);
    my $i1 = at($idx, $n - 2 - $j + 1);
    $p($i) .= $p($i1) if $p($i1) < $p($i);
  }
  
  return $p
}

##################################################

=item BenjaminiHochbergLimit

  my ($plim, $ns) = BenjaminiHochbergLimit($p, $alpha);
  
For given significance $alpha, calculates limit for Benjamini-Hochberg ($plim) correction and the number of significant elements ($ns). If there are no significant p-values ($ns = 0), Bonferroni limit ($alpha/$n) is returned.

This subroutine is B<much> faster than BenjaminiHochberg (see above).

=cut 

sub BenjaminiHochbergLimit
{
  my ($p_, $alpha) = @_;
  
  $alpha ||= 0.05;
  
  my $p = $p_->copy();
  my $n = nelem($p);
  my $ps = qsort $p;
  my $k = 1;
  while($k <= $n && at($ps, $k-1) * $n / $k <= $alpha)
    {$k++}
  #$k--;
  #my $plim = ($k >= 0) ? at($ps, $k) : $alpha/$n;
  my $plim = ($k >= 1) ? $k / $n * $alpha : $alpha/$n;
  return $plim unless wantarray;
  return $plim, $k
}

##################################################

=item HolmBonferroni

  my $pc = HolmBonferroni($p);

Corrects a list of p-values (piddle $p) for multiple testing using step-down Holm-Bonferroni correction. Returns a piddle with corrected p-values.

=cut

sub HolmBonferroni
{
  my $p_ = shift;
  my $p = $p_->copy();
  
  my $n = nelem($p);
  my $idx = qsorti $p;
  my $c = $n - sequence($n);
  $p($idx) *= $c;
  for my $j (1 .. $n - 1)
  {
    my $i = at($idx, $j);
    my $i1 = at($idx, $j - 1);
    $p($i) .= $p($i1) if $p($i1) > $p($i);
  }
  
  $p->where($p>1) .= 1;
  return $p
}

##################################################

=item HolmBonferroniLimit

  my $plim = HolmBonferroniLimit($p, $alpha);
  my ($plim, $ns) = HolmBonferroniLimit($p, $alpha);
  
For given significance $alpha, calculates limit for Holm-Bonferroni ($plim) correction and the number of significant elements ($ns). If there are no significant p-values ($ns = 0), Bonferroni limit ($alpha/$n) is returned.

This subroutine is B<much> faster than HolmBonferroni (see above).

=cut 

sub HolmBonferroniLimit
{
  my ($p_, $alpha) = @_;
  
  $alpha ||= 0.05;
  
  my $p = $p_->copy();
  my $n = nelem($p);
  my $ps = qsort $p;
  my $k = 0;
  while($k < $n && at($ps, $k) * ($n - $k + 2) <= $alpha)
    {$k++}
  $k--;
  my $plim = ($k >= 0) ? at($ps, $k) : $alpha/$n;
  
  return $plim unless wantarray;
  return $plim, $k + 1
}

sub HolmBonferroniLimit_
# Alas, this doesn't work with data with repeated p-values
# because vsearch needs monotonic increase
{
  my ($p_, $alpha) = @_;
  
  $alpha ||= 0.05;
  
  my $p = $p_->copy();
  my $n = nelem($p);
  my $ps = qsort $p;
  my $k = vsearch($alpha, $ps * ($n - sequence($n))) - 1;
  print $ps*($n-sequence($n)), "\n";
  
  my $plim = ($k >= 0) ? at($ps, $k) : $alpha/$n;
  return $plim, $k + 1
}

##################################################

=item PoissonTest

  my $p = PoissonTest($x, $type);

Performs test for Poisson distribution. There are two tests avaialable:

  $type = 'dispersion' (default) - dispersion test
  $type = 'logratio' - log-ratio test

=cut

sub PoissonTest
{
  my ($x, $type) = @_;
  
  $type ||= 'dispersion';
  
  my $n = nelem($x);
  my ($m, $s) = stats($x);
  my $chi2;
  if($type eq 'logratio')
    {$chi2 = 2 * sum($x * log($x/$m))}
  elsif($type eq 'dispersion')
    {$chi2 = sum(($x - $m)**2 / $m)}
  else
    {die "Unrecognized Poisson test type $type\n"}
    
  my $p = 1 - pchisq($chi2, $n - 1);
  
  return $p
}

##################################################

=item NBKSTest

  my ($D, $p) = NBKSTest($x, $sim, $R);
  
Kolmogorov-Smirnov goodness-of-fit test for negative binomial distribution. Uses discrete version of the KS test, based on Conver (1972) and Glesser (1985). The test is done by R package developed by Arnold & Emerson (2011).

  $x - data piddle
  $R - (optional) path to R, where 'dgof' is installed
  
Conover WJ (1972), "A Kolmogorov goodness-of-fit test for discontinuous distributions", J Am Stat Assoc, 67, 591

Gleser LJ (1985), "Exact power of goodness-of-fit tests of Kolmogorov type for discontinuous distributions", J Am Stat Assoc, 80, 954

Arnold TB, Emerson JW (2011), "Nonparametric goodness-of-fit tests for discrete null distributions", The R journal, 3, 34

=cut

sub NBKSTest
{
  my ($x, $sim, $R) = @_;
  $R ||= 'R';
  
  my ($loglik, %pars) = FitDistr(floor($x+0.5), 'negative binomial');
  my $mu = $pars{mu}[0];
  my $size = $pars{size}[0];
  
  my $max = ceil(max($x) * 2);  # good enough?
  
  #my $dir = '/tmp';
  my $dir = (defined $ENV{TMPDIR}) ? $ENV{TMPDIR} : '/tmp';
  
  
  my ($FR, $Rfile) = tempfile('KSNBXXXXXX', DIR=>$dir, SUFFIX=>'.R', UNLINK=>1); 
  #my ($Fout, $outfile) = tempfile('KSNBXXXXX', DIR=>$dir, SUFFIX=>'.out');
  #close $Fout; 
  
  my $outfile = mktemp("$dir/KSNBXXXXXX");
  
  #unlink $infile, $outfile, $Rfile;
  
  my $cx = join ',', list($x);
  my $sm = ($sim) ? ", simulate.p.value=T, B=$sim" : '';

  print $FR <<EOF;
library("dgof")
x <- c($cx)
f <- stepfun(0:$max, c(0, pnbinom(0:$max, mu=$mu, size=$size))) 
dgof::ks.test(x, f$sm)
EOF
  close $FR;
  call("$R --vanilla < $Rfile >& $outfile", my $stat, 1);
  #die "Problem with R call fitdistr\n"if !-e $outfile || $stat;
  if(!-e $outfile || $stat)
  {
    warn "Problem with KSNBTest\n";
    return (0, ());
  }
  
  local *F;
  open F, $outfile or die "Cannot open $outfile\n";
  my ($D, $p);
  while(my $line = <F>)
  {
    chomp $line;
    if($line =~ /D\s+=\s+(.*)\,\s+p-value\s+=\s+(.*)/)
    {
      $D = $1;
      $p = $2;
      last
    }
  }
  die "Cannot parse R output\n" unless defined $D && defined $p;
  close F;

  unlink $Rfile, $outfile;
  return $D, $p;
}

##################################################

=item NBBootstrapTest

  my ($p, $S) = NBBootstrapTest($stat, $x, $ncrit, $maxn);

Test for negative binomial distribution, using either Meintanis or Anderson-Darling statistic. 

  $stat - either 'm' or 'a'
  $x - piddle with data to test
  $ncrit - required minimum number of crtical tests (default 30)
  $maxn - maximum number of simulations (default 1e7)

=cut

sub NBBootstrapTest
{
  my ($stat, $x, $ncrit, $maxn, $verbouse) = @_;
  
  $ncrit ||= 30;
  $maxn ||= 1e7;
  
  my $n = nelem($x);  
  my ($m, $s) = stats($x);
  my $S = ($stat eq 'm') ? 
    MeintanisStatistic($x) : 
    AndersonDarlingNBStatistic($x);

  my ($loglik, %pars) = FitDistr(floor($x+0.5), 'negative binomial');
  ($m, $s) = GetRMS('negative binomial', %pars);
  #print "$m $s    $M $S  ";
  return 0 if $s**2 <= $m;
  my ($bx, $by) = CreateNBinomialTables($m, $s);
  my $crit = 0;
  my $tot = 1;
  
  while($crit < $ncrit)
  {
    printf "\b\b\b\b\b\b\b\b\b\b\b\b%6d %4d", $tot, $crit if $verbouse && $tot % 1000 == 0;
    my $xsim = GenerateNBinomial($bx, $by, $n);
    my $Ssim = ($stat eq 'm') ? 
      MeintanisStatistic($xsim) : 
      AndersonDarlingNBStatistic($xsim, $m, $s); 
    $crit++ if $Ssim > $S;
    $tot++;
    last if $tot > $maxn;
  }
  print "\n" if $verbouse;
  my $p = $crit / $tot;
  #print "N=$tot  ";

  return $p;
}

##################################################

=item MeintanisStatistic

  my $T = MeintanisStatistic($x, $a);

Statistic used for negative binomial distribution test, see Meintanis (2005), "Transform method for testing the negative binomial hypothesis", Statistica, 3, 294.

  $x - piddle with data to test
  $a - a parameter, default value is 5

=cut

sub MeintanisStatistic
{
  my ($x, $a) = @_;
  $a ||= 5;
  
  my $n = nelem($x);
  my ($m, $s1, $med, $min, $max, $adev, $s) = stats($x);
  my $rho = ($s**2 - $m) / $m;
  my $xp = $x + transpose($x);
  
  #print "x=$x\n";
  #print "xp=$xp\n";
  #print "I(xp+a) = ", I($xp+$a), "\n";
  #print "sum=", sum(I($xp+$a)), "\n";
  
  my $t1 = $m**2 * sum(I($xp + $a));
  #print "t1=$t1\n";
  
  my $t21 = (1 + $rho) * sumover(I($xp + $a - 1));
  my $t22 = -$rho * sumover(I($xp + $a));
  my $t2 = -2*$m*sum($x * ($t21 + $t22));
  #print "t2=$t2\n";
  
  my $t31 = (1 + $rho)**2 * I($xp + $a - 2);
  my $t32 = $rho**2 * I($xp + $a);
  my $t33 = -2 * $rho * (1 + $rho) * I($xp + $a - 1);
  my $t3 = sum($x * sumover($x * ($t31 + $t32 + $t33)));
  #print "t3=$t3\n";
  my $T = ($t1 + $t2 + $t3) / $n;
}

sub I
{
  my $b = shift;
  1 / (1 + $b)
}

##################################################

=item AndersonDarlingNBStatistic

  my $A2 = AndersonDarlingNBStatistic($x, $M, $S);

Anderson-Darling statistic, see Anderson, T.W., and D.A. Darling (1952), “Asymptotic theory of certain ‘goodness-of-ﬁt’ criteria based on stochastic processes,” Ann. Math, Statistics, Vol. 23, 193–212.

This one is calculated for the negative binomial distribution. Parameters (mean $M and standard deviation $S) of the NB distribution can be given or, if not present, will be estimated from the maximum likelihood fit. The distribution of this statistic is not known in closed form, so the critical values/probabilities have to be calculated from simulation. See NBBootstrapTest

=cut

sub AndersonDarlingNBStatistic
{
  my ($x, $M, $S) = @_;
  
  unless(defined $M && defined $S)
  {
    my ($loglik, %pars) = FitDistr(floor($x+0.5), 'negative binomial');
    ($M, $S) = GetRMS('negative binomial', %pars);
  }
  
  my $xs = qsort $x;
  my $n = nelem($xs);
  my $s = 0;
  for my $k (1 .. $n)
  {
    my $F1 = nbinomms(at($xs, $k-1), $M, $S);
    my $F2 = nbinomms(at($xs, $n-$k), $M, $S);
    $s += (2*$k - 1)/$n * (log($F1) + log(1 - $F2)); 
  }
  return -$n - $s
}

##################################################

=item HinzGurland

  my ($Xf, $p) = HinzGurland($x);

Test for negative binomial distribution. Provides with the chi-square statistics, $Xf and p-value, $p.

This is implementation of the algorithm from Hinz & Gurland (1970), "A test of fit for the negative binomial and other contagious distributions", J. Am. Stat. Assoc., 65, 887.

Imiportant note: quantities in this test are calculated based on the zero frequency (p0), i.e., normalized number of zeroes in the sample. If p0 = 0 (there are no zeroes in the sample), the Jacobian matrix (A5 in the paper) doesn't seem to make much sense. Matrix \Sigma is singular and cannot be inversed, so the test fails.

=cut


sub HinzGurland
{
  my ($x, $devel) = @_;
  
  my $n = nelem($x);
  print "\n\nData=$x\n" if $devel;
  
  return (0, 1) if max($x) < 3;
  
  my $m1 = sum($x) / $n;
  my $m2 = sum($x**2) / $n;
  my $m3 = sum($x**3) / $n;
  my $m4 = sum($x**4) / $n;
  
  my $k1 = $m1;
  my $k2 = $m2 - $m1*$k1;
  my $k3 = $m3 - $m2*$k1 - 2 * $m1*$k2;
  my $k4 = $m4 - $m3*$k1 - 3 * $m2*$k2 - 3*$m1*$k3;
  
  my $f1 = $k1;
  my $f2 = $k2 - $k1;
  my $f3 = $k3 - 3*$k2 + 2*$k1;
  my $f4 = $k4 - 6*$k3 + 11*$k2 - 6*$k1;

  if($devel)
  {
    print "m1=$m1\nm2=$m2\nm3=$m3\nm4=$m4\n";
    print "k1=$k1\nk2=$k2\nk3=$k3\nk4=$k4\n";
    print "f1=$f1\nf2=$f2\nf3=$f3\nf4=$f4\n";
  }
  
  my $h0 = $f1;
  my $h1 = $f2 / $f1;

  return 0, 0 if $h1 == 0;
  
  #my ($loglik, %pars) = FitDistr($x, 'negative binomial');
  #my ($M, $S) = GetRMS('negative binomial', %pars);
  my ($M, $S) = stats($x);
  my $q = $S*$S - $M;
  my $rf = $M*$M / $q;
  my $pf = $q / ($S*$S);    # best fitting p and r
  my $p0 = (1 - $pf)**$rf;
  #my $p0 = nelem($x->where($x==0)) / $n;  # this is not very good
  #print "$p0  $p00\n";  
  
  my $A = (1 + $h1)**(($h0 + $h1)/$h1);
  my $a1 = $p0 * $A - 1;
  my $t = pdl($h0, $h1, $a1)->transpose();
  my $Wp = pdl([1, 0, 0], [0, 1, 1]);
  my $W = transpose($Wp);

  print "h0=$h0\nh1=$h1\np0=$p0\na1=$a1\nt=$t\nWp=$Wp\n" if $devel;
  
  my $xi1 = $k1*$k3 - $k2**2;
  my $xi2 = $k1*($k1*($k4 + 2*$k2**2) - $k2*$k3) - $k2*($k1*$k3 - $k2**2);
  my $Omega = pdl(
    [$k2,         $xi1/$k1**2, -$p0*$k1],
    [$xi1/$k1**2, $xi2/$k1**4,  $p0*$k1],
    [-$p0*$k1,    $p0*$k1,      $p0*(1 - $p0)]
  );
  
  print "xi1=$xi1\nxi2=$xi2\nOmega=$Omega\n" if $devel;
  
  my $j13 = ($a1 + 1) * log(1 + $h1) / $h1;
  my $j23 = ($a1 + 1) * (($h0 + $h1)/($h1*(1 + $h1)) - $h0*log(1 + $h1)/$h1**2);
  #my $j33 = ($a1 + 1) / $p0;
  my $j33 = $A; 
  my $J = pdl(
    [1,       0,    0],
    [0,       1,    0],
    [$j13, $j23, $j33]
  );
  my $Jp = transpose($J);
  my $Sigma = $J x $Omega x $Jp;
  #$Sigma(2,2) .= 1e-32 if $p0 == 0;
  
  print "Jp=$Jp\nSigma=$Sigma\n" if $devel;

  my $Sigma1 = inv($Sigma, {s=>1});
  #print "Sigma singular\n" unless defined $Sigma1;
  return 0, 0 unless defined $Sigma1;
  
  print "Sigma1=$Sigma1" if $devel;
  
  my $R = $Wp x $Sigma1 x $W;
  my $Theta = inv($R) x $Wp x $Sigma1 x $t;
  return 0, 0 unless defined $Theta;

  print "Theta=$Theta\n" if $devel;

  my $Q = $t - $W x $Theta;
  my $Qp = transpose($Q);
  my $Xf_ = flat($n * $Qp x $Sigma1 x $Q);
  my $Xf = at($Xf_, 0);
  #print "Xf < 0\n" if $Xf < 0;
  return 0, 0 if $Xf < 0;
  
  my $p = 1 - pchisq($Xf, 1);
  
  print "Xf=$Xf\np=$p\n" if $devel;
  return $Xf, $p;
}


##################################################

=item NBMiTest

  my $p = NBMiTest($x, %opt);

Test for negative binomial distribution, Mi G., Di Y., Schafer D.W., PLOS One (2015) 

  $x - piddle with data to test

  $opt{model} - model to use, default is NB2
  $opt{nsim} - number of MC simulations, default is 999
  $opt{method} -  ML or APL

=cut

sub NBMiTest
{
  my ($x, %opt) = @_;
  
  my %defaults = (
    model => 'NB2',
    nsim => 999,
    method => 'ML',
    fname => undef,
    debug => undef
  ); 
  %opt = parse(\%defaults, \%opt);
  my ($model, $nsim, $method, $fname, $debug) = map {$opt{$_}} qw(model nsim method fname debug);
  
  my $n = nelem($x);
  my $c = 'c(' . join(',', list($x)) . ')';
  
  my $tmpdir = './mitmp';
  mkdir $tmpdir, 0777 unless -d $tmpdir;
  
  my $f = (defined $fname) ? "_${fname}_" : '';  # for debugging
  
  my $outfile = "$tmpdir/" . mktemp("result${f}XXXXXX");
  my $Rfile = "$tmpdir/" . mktemp("Rscript${f}XXXXXX");
  my $errfile = "$tmpdir/" . mktemp("error${f}XXXXXX");
  unlink $outfile, $Rfile;

  local *R;
  open R, ">$Rfile" or die "Cannot open temporary file $Rfile\n";
  print R <<EOF;
library(NBGOF)
res = nb.gof.v($c, matrix(rep(1,$n)), sim=$nsim, model="$model", method="$method")
res\$new.pval
EOF
  close R;
  #call("R --vanilla < $Rfile > $outfile", my $stat);
  #`R --vanilla < $Rfile > $outfile 2> /dev/null`;
  `R --vanilla < $Rfile > $outfile 2> $errfile`;
  die "Problem with R call nls\n"if !-e $outfile;

  local *F;
  open F, $outfile or die;
  my $line;
  while($line = <F>)
    {last if $line =~ /^\[1\]/}
  chomp $line;
  my ($ss, $p) = split(" ", $line);
  
  unlink $outfile, $Rfile unless $debug;
  
  return $p;
}

