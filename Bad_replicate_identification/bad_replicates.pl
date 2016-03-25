#!/sw/bin/perl

=head1 NAME

B<bad_replicates.pl>

=head1 DESCRIPTION

Creates a combined figure to find bad replicates. Adds bad replicate numbers, as defined in C<defs.dat> file.

Requires pileup statistics created with C<pileup_var.pl>.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use GetOpt;
use GRNASeq;

use Stats;
use PLgraphs;
use Tools;

#require 'mylib.pl';

Stamp();
ReadDefs();

$| = 1;
my $PI = 4 * atan(1);

my $limit = 5;
my $ntrim = 3;
my $minm = 0;
my $method = 's';   # c - Chauvenet, s - sigma
my $off = 0.07;
GetMyOptions(
  'limit=f' => \$limit,
  'ntrim=i' => \$ntrim,
  #'minm=f' => \$minm,      # not used
  #'method=s' => \$method,  # not used
  'off=f' => \$off,
);


my @cl = split(/:/, $CleanExclude);
my %clean = ();
for my $i (0 .. 1)
{
  my @c = split(/,/, $cl[$i]);
  $clean{$conds[$i]}{$_} = 1 for (@c);
}


my %dat = GetNormalizedData(
  nonzero=>1,
  norm=>$norm,
  type=>$type
);


my $w = NewWindow($psfile);
my $cs = 0.6;
    
$PL_ymin = 0.08;
$PL_ymax = 0.95;
$PL_xmax = 0.95;
$PL_deltax = 0.02;
$PL_deltay = 0.02;
$PL_nx = 2;
$PL_ny = 4;


my %opt = ();
my $px = 1;
for my $cond (@conds)
{
  my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel) = @{$dat{$cond}};
  print "$ngen genes and $nrep replicates read.\n";
  #my %g2i = GenHash($genes, $sel);
  
  ### median correlation
  
  my $C = CorrelationMatrix($d);
  my ($MC, $MS, $r) = statsover($C);
  print statsover($C);
  PlotReps($w, $px, 1, $r, $cond, top=>1, min=>0.65, max=>1, ylab=>'Median correlation');
  TopText($w, $cond, y=>2, x=>0);
  
  ### outliers
  
  my $fout = CountOutliers($d);
  PlotReps($w, $px, 2, $fout, $cond, min=>1e-4, max=>0.45, ylab=>'Outlier fraction');

  ### Median chi-squared for non-uniformity of read distribution
  
  my ($chisq, $se) = PileupChisq($cond, $nrep);
  PlotReps($w, $px, 3, $chisq, $cond, min=>0, max=>30, ylab=>'Median chisq');
  
  ### SCORE
  
  my $score = log10((1 - $r) * $fout * ($chisq - 1));
  PlotReps($w, $px, 4, $score, $cond, min=>-4, max=>1, y0=>-5, ylab=>'Score');

  #TopText($w, $cond, x=>0, y=>1);
  
  $px++ 
}
$w->close();


##################################################

sub PlotReps
{
  my ($w, $px, $py, $y, $cond, %opt) = @_;
  
  my ($min, $max) = ($opt{min}, $opt{max});
  ($min, $max) = BoxMinMax($y) unless defined $min && defined $max;
  my $offset = $opt{offset};
  $offset = $off * ($max - $min) unless defined $offset;
  $offset *= -1 if $opt{top};
  
  my $nrep = nelem($y);
  PlotPanelXY($w, $px, $py,
    BOX=>[0, $nrep+1, $min, $max],
    xbox=>'I', ybox=>'I',
    XLAB => 'Replicate', YLAB => $opt{ylab},
    CHARSIZE => $cs
  );
  
  my $y0 = ($opt{top}) ? 1 : 0;
  $y0 = $opt{y0} if defined $opt{y0};
  PointPlot($w, sequence($nrep)+1, $y, SYMBOLSIZE=>1);
  for my $i (1 .. $nrep)
  {
    my $yy = at($y, $i - 1);
    LinePlot($w, pdl($i,$i), pdl($y0, $yy), COLOR=>'GREY');
    if($clean{$cond}{$i})
      {$w->text($i, CHARSIZE=>0.4, COLOR=>'RED', TEXTPOSITION=>[$i, $yy+$offset, 0, 0, 0.5])}
  }
  
  if($px == 1)
    {TopText($w, '(' . chr(96+$py) . ')', x=>0.02, CHARSIZE=>$cs)}
}

##################################################

sub CorrelationMatrix
{
  my ($d) = @_;
  
  my ($nrep, $ngen) = dims($d);
  
  my $C = zeroes($nrep, $nrep) + 1;
  for my $i (0 .. $nrep - 2)
  {
    for my $j ($i + 1 .. $nrep - 1)
    {
      my $r1 = $d($i,;-);
      my $r2 = $d($j,;-);
      my $sel = which(($r1>0) & ($r2>0));
      my $r = corr($r1($sel), $r2($sel));
      #printf "%02d %02d %8.6f\n", $i+1, $j+1, $r;
      $C($i, $j) .= $r;
      $C($j, $i) .= $r;
    }
  }
  return $C;
}

##################################################

sub CountOutliers
{
  my ($d) = @_;
  
  my ($nrep, $ngen) = dims($d);
  my $up = zeroes($nrep);
  my $down = zeroes($nrep);
  my %out = ();
  print "Processing outliers...       ";
  for my $g (0 .. $ngen - 1)
  {
    Percent($g/($ngen-1));
    my $dat = $d(,$g)->flat();
    #my $md = at($m, $g);
    #next if $md < $minm;

    my ($clean_dat, $iup, $idown, $M, $S) = RejectOutliersSigma($dat, $limit, $ntrim);
    $up($iup)++;
    $down($idown)++;
  }
  print "\b\b\b\b\b\b\b done   \n";
  
  $up /= $ngen;
  $down /= $ngen;
  my $all = $up + $down;
}

##################################################

sub PileupChisq
{
  my ($cond, $nrep) = @_;
    
  my $mfile = "$pileupdir/${cond}_pileupvar_mean.tsv";
  my $sfile = "$pileupdir/${cond}_pileupvar_sd.tsv";
  
  #local *F;
  #open F, $mfile or die "Cannot open $mfile\n";
  #my $line = <F>;
  #my @s = split " ", $line;
  #close F;
  #my $nrep = @s - 1;
      
  print "Reading pileup files for $cond...";
  my @M = rcols $mfile, (1 .. $nrep); 
  my @S = rcols $sfile, (1 .. $nrep);
  print " done\n";

  my @repvar = ();
  my @se = ();
  for my $i (0 .. @M - 1)
  {
    my $M = $M[$i];
    my $S = $S[$i];
    my $sel = which($M > 0);

    my $CV = $S($sel) / $M($sel);
    my $CQ = $S($sel)**2 / $M($sel);
   
    my ($m, $s, $med) = stats($CQ);
    my $qCQ = qsort $CQ;
    my $perc = at($qCQ, 0.9*nelem($qCQ));
    
    my $n = nelem($CQ);
    my $sCQ = qsort $CQ;
    my $L = floor($n/2) - ceil(sqrt($n/4));
    my $U = $n - $L;
    
    # standard error of the median
    my $se = (at($sCQ, $U-1) - at($sCQ, $L)) / 2;
    push @se, $se;
    push @repvar, $med;
    #push @repvar, $perc;
  }  
  my $repvar = pdl(@repvar);
  my $se = pdl(@se);
  
  return $repvar, $se;
}



=head1 SYNOPSIS

  bad_replicates.pl  
    
=head1 OPTIONS

=over 4

=item B<-limit>=I<number>

Significance limit for outlier selection. This is a limit of number of standard deviations from the trimmed mean. The default value is 5.

=item B<-ntrim>=I<number>

Number of data points to remove on each side for the trimmed mean. The default value is 3.

=item B<-off>=I<number>

Offset for plotting bad replicate numbers above/below data points.

=item B<-psfile>=I<name>

Name of the postsript file to redirect output to.

=back
