#!/sw/bin/perl

=head1 NAME

B<plot_ms>

=head1 DESCRIPTION

Plot mean, standard deviation and a few over things. This script can make lots of plots that we used to investigate data.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use GRNASeq;
use GetOpt;

use Stats;
use PLgraphs;
use Tools;
use Distribution;
use Switch 'Perl6';

$| = 1;
Stamp();
ReadDefs();

my $graph = 'mv2';
my ($set1, $set2, $withtest, $withfull, $withmean, $withloess, $sellist);
my $dispersion = -1;
my $testtype = 'dispersion';
my $plotms;
my $outrep;
my $nvar;
my $log;
my $witherr;
my $cs = 0.7;
GetMyOptions(
  #'graph=s' => \$graph,
  #'set1=s' => \$set1,
  #'set2=s' => \$set2,
  #'testtype=s' => \$testtype,
  #'outrep=i' => \$outrep,
  #'nvar=i' => \$nvar,
  #withtest => \$withtest,
  #withfull => \$withfull,
  #witherr => \$witherr,
  #withmean => \$withmean,
  withloess => \$withloess,
  #plotms => \$plotms,
  #'dispersion=f' => \$dispersion,
  #log => \$log,
  #'cs=f' => \$cs,
  #'sellist=s' => \$sellist,
);

my ($minm, $maxm) = (-2, 6);
$PL_xmax = 0.9;
$PL_ymax = 0.9;


my %dat = GetNormalizedData(
  replicate=>$rep,
  exclude=>$exclude,
  include=>$include,
  clip=>$clip,
  nonzero=>$nonzero,
  type=>$type,
  norm=>$norm,
  clean=>$clean,
  spclean=>$spclean,
  randrep=>$randrep,
);

my @hl = ();                       # highlight
my %g2i;
my $N;
my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel, $f, $xsel);
my ($lm, $ls, $full);
unless($graph eq 'mv2' || $graph eq 'comp')
{
  ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel, $f, $xsel) = @{$dat{$cond}};
  %g2i = GenHash($genes, $sel);
  $N = sumover($d>0);
  $lm = log10($m);
  $ls = log10($s);
  print " done\n";
  print "$ngen non-zero ($nonzero) genes in $nrep replicates read.\n";
  $full = which($N == $nrep);
  @hl = sellist(\%g2i);
}

if(defined $outrep)
{
  $genlist = OutliersFile(undef, $cond, $outrep, $norm);
  $name .= " with outliers rep$outrep";
}

my $glist = GeneIdList($genlist, \%g2i) if defined $genlist;

#print "$glist\n";
#print "f=$f\n";

given($graph)
{
  when 'ms'   {PlotMeanSD(0, $withtest)}
  when 'mv'   {PlotMeanSD(1, $withtest)}
  when 'mv2' {PlotMeanVar2()}
  when 'mvm'   {PlotMeanVM()}
  when 'mlsm' {PlotMeanlSM()}
  when 'mo'   {PlotMeanOver()}  
  when 'cv'   {CompareCumulVar($set1, $set2)}
  when 'vs' {PlotVarshift(1)}
  when 've' {PlotVarErr()}
  when 'comp' {Composition()}
  when 'do' {PlotDistOver()}
}


##############################################

sub PlotMeanVar2
{
  my $win = NewWindow($psfile);

  $PL_nx = @conds;
  $PL_ymin = 0.5;
  $PL_deltax = 0.04;
  my $pan = 1;
  for my $cond (@conds)
  {
    ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel) = @{$dat{$cond}};
    my %g2i = GenHash($genes, $sel);
    my @hl = sellist(\%g2i);
    
    my $o = ($s**2 - $m) / $m**2;
    my $x = log10($m);
    my $y = 2 * log10($s);
    
    my $dm = $s / sqrt($nrep);                   # SE of the mean
    my $dv = $s**2 * sqrt(2 / ($nrep - 1));  # SE of the variance
    
    #print min($y), " ", max($y), "\n";
    
    my $pcol = ($sellist) ? 'GREY' : 'BLACK';
    
    PlotPanelN($win, $pan++,
      BOX => [$minm, $maxm, -2.3, 10],
      XLAB => 'log mean',
      YLAB => 'log variance',
      CHARSIZE => $cs,
    );
    
    if($witherr)
    {
      my $x1 = log10($m - $dm);
      $x1->where($m - $dm < 1e-3) .= -3;
      my $x2 = log10($m + $dm);
      my $x0 = ($x1 + $x2) / 2;
      my $dx = $x2 - $x1;
      my $y1 = log10($s**2 - $dv);
      my $y2 = log10($s**2 + $dv);
      my $y0 = ($y1 + $y2) / 2;
      my $dy = $y2 - $y1;
      
      #my $sel = which($x < 0);
      #wcols $m($sel), $dm($sel), $x0($sel), $x1($sel), $x2($sel);
      
      PointPlot($win, $x0, $y, SYMBOLSIZE=>0.001, MINTICKSIZE=>0.001, COLOR=>'GREY', XERRORBAR=>$dx);
      PointPlot($win, $x, $y0, SYMBOLSIZE=>0.001, MINTICKSIZE=>0.001, COLOR=>'GREY', YERRORBAR=>$dy);
    }
    
    PointPlot($win, $x, $y, SYMBOLSIZE=>0.3, COLOR=>$pcol);
    TopText($win, $cond);
    #print min($y), " ", max($y), "\n";
    
    my $xx = pdl($minm, $maxm);
    my $yy = $xx;
    LinePlot($win, $xx, $yy, COLOR=>'RED', LINEWIDTH=>2, LINESTYLE=>4);  # poisson
    LinePlot($win, pdl(-10,-10), pdl(-10,-10), LINESTYLE=>1);
    PlotRunningMean($win, $m, $s);
    PlotDispersionLine($win, $m, $s);
    PlotLoess($win, $x, $y);
    
    for my $i (0 .. @hl - 1)
      {PointPlot($win, $x($hl[$i]), $y($hl[$i]), SYMBOLSIZE=>1.3, COLOR=>$colours[$i])}
  
  #for my $phi (0.001,  0.1)
  #  {displine($win, $phi, 2, [215,215,255], 3)}  
  LinePlot($win, pdl(-10,-10), pdl(-10,-10), LINESTYLE=>1);
  }

  $win->close();  
}


sub PlotDispersionLine
{
  my ($win, $m, $s, $f) = @_;
  $f ||= 2;  
  
  return unless $dispersion;
  my $o = ($s**2 - $m) / $m**2;
  my $disp = $dispersion;
  if($dispersion < 0)
  {
    #$disp = median($o);
    my ($ml, $mu);
    ($disp, $ml, $mu) = MedianCI($o, 0.95);
    $ml -= $disp; $mu -= $disp;
    print "Median dispersion = $disp [$ml $mu]\n";
       
  }
  displine($win, $disp, $f);
}

sub displine
{
  my ($win, $disp, $f, $col, $styl) = @_;
  $f ||= 2;
  
  $col ||= 'GOLD2';
  $styl ||= 1;
  
  my $lmm = ($maxm - $minm) * sequence(100)/100 + $minm;
  my $mm = 10**$lmm;
  my $so = $f * 0.5 * log10($mm * ($disp * $mm + 1));
  LinePlot($win, $lmm, $so, COLOR=>$col, LINEWIDTH=>2, LINESTYLE=>$styl);
}

sub PlotRunningMean
{
  my ($win, $m_, $s_, $winsize, $o) = @_;
  $winsize ||= 30;
  
  my $idx = qsorti $m_;
  my $m = $m_($idx);
  my $s = $s_($idx);
  
  return unless $withmean;
  my $n = nelem($m);
  my $v = ($o) ? $s : $s**2;
  my @x = (); my @y = ();
  for my $i (0 .. $n - 1)
  {
    my $i1 = $i - $winsize/2;
    my $i2 = $i + $winsize/2;
    if($i1 < 0)
    {
      $i2 = $i2 + $i1;
      $i1 = 0
    }
    if($i2 > $n - 1)
    {
      $i1 = $i1 + ($i2 - $n + 1);
      $i2 = $n - 1
    }
    #my $x = mean(log10($m($i1:$i2)));
    my $x = ($o) ? at($m, $i) : log10(at($m, $i));
    my $y = ($o) ? mean($v($i1:$i2)) : mean(log10($v($i1:$i2)));
    push @x, $x;
    push @y, $y;
    $i1++;
    $i2++;
  }
  
  LinePlot($win, pdl(@x), pdl(@y), COLOR=>[180,255,180], LINEWIDTH=>2);
}

##############################################


sub PlotMeanSD
{
  my ($var, $test) = @_;
  
  my $win = NewWindow($psfile);
  
  my ($miny, $maxy) = ($var) ? (-1.9, 9) : (-1, 5);
  my $ylab = ($var) ? 'Variance' : 'SD';
  my $f = ($var) ? 2 : 1;
  #my $ybox = ($var) ? 'L' : 'L';
  
  my $x = $lm;
  my $y = $f * $ls;
  
  my %box = (xbox => 'L', ybox => 'L') unless $log;
  
  PlotPanelN($win, 1,
    BOX => [$minm, $maxm, $miny, $maxy],
    XLAB => 'Mean',
    YLAB => $ylab,
    CHARSIZE => $cs,
    %box
  );
  
  my $col = ($test || $withfull || $genlist || $sellist) ? 'GREY' : 'BLACK';
  PointPlot($win, $x, $y, SYMBOLSIZE=>0.3, COLOR=>$col);
  if($test)
  {
    my $P = PoissonFromData($d, $testtype);
    my $Plim = HolmBonferroniLimit($P);
    my $sel = which($P < $Plim);
    
    PointPlot($win, $x($sel), $y($sel), SYMBOLSIZE=>0.3, COLOR=>'RED');
  }
  PlotDispersionLine($win, $m, $s, $f);
  
  if($withfull)
    {PointPlot($win, $x($full), $y($full), SYMBOLSIZE=>0.3, COLOR=>'BLACK')}
  
  if($genlist)
    {PointPlot($win, $x($glist), $y($glist), SYMBOLSIZE=>0.3, COLOR=>'BLACK')}
    
  for my $i (0 .. @hl - 1)
    {PointPlot($win, $x($hl[$i]), $y($hl[$i]), SYMBOLSIZE=>1.3, COLOR=>$colours[$i])}
  

  my $xx = 9*sequence(1000)/1000 - 2;
  my $yy = 0.5 * $f * $xx;
  LinePlot($win, $xx, $yy, COLOR=>'RED', LINEWIDTH=>2);  # poisson
  #LinePlot($win, $xx, $f*$xx, COLOR=>'GOLD2', LINEWIDTH=>1);
  
  #while(1)
  #{
  #  my %gin = plGetCursor();
  #  my $x = $minm + ($maxm-$minm)*($gin{dX} - $PL_xmin) / ($PL_xmax - $PL_xmin);
  #  my $y = $miny + ($maxy-$miny)*($gin{dY} - $PL_ymin) / ($PL_ymax - $PL_ymin);
  #  printf "%.3g %.3g\n", 10**$x, 10**$y;
  #}
  if($plotms)
  {
    my ($x, $y) = rcols 'ms.dat';
    wcols $x, $y;
    PointPlot($win, log10($x), log10($y), COLOR=>'RED', SYMBOL=>851, SYMBOLSIZE=>0.8);
  }
  
  
  Title($win);
  $win->close();
}

##############################################

sub PlotVarshift
{
  my $win = NewWindow($psfile);
  
  $PL_nx = 2;
  $PL_ny = 2;
  $PL_deltax = 0.07;

  my $x = $lm;
  $nvar ||= 3;
  my ($merr, $verr) = Varsim($nvar);

  PlotPanelN($win, 1,
    BOX => [$minm, $maxm, 0, max($merr) * 1.1],
    XLAB => 'log Mean',
    YLAB => 'Fractional error in mean',
    XBOX => 'BNSTI', YBOX => 'BNSTI',
    CHARSIZE => $cs,
    forcexlab=>1, forceylab=>1
  );

  PointPlot($win, $x, $merr);

  PlotPanelN($win, 2,
    BOX => [$minm, $maxm, 0, max($verr) * 1.1],
    XLAB => 'log Mean',
    YLAB => 'Fractional error in variance',
    XBOX => 'BNSTI', YBOX => 'BNSTI',
    forcexlab=>1, forceylab=>1,
    CHARSIZE => $cs,
  );

  PointPlot($win, $x, $verr);
  
  my $f = sqrt(2 / ($nvar - 1));
  LinePlot($win, pdl($minm, $maxm), pdl($f, $f), COLOR=>'RED');
  
  FullBox($win);
  TopText($win, "Mean and variance error for $nvar replicates", c=>0.5, x=>0.5, y=>-2);
  TopText($win, $name, c=>0.5, x=>0.5, y=>-4);
  
  $win->close();
}

sub PlotVarErr
#
# This plot doesn't show any information about data. It only shows
# that the fractional error of the variance is (always) sqrt(2/(n-1)).
# This is trivial.
#
{
  my $win = NewWindow($psfile);
  
  my $maxn = $nrep - 1;
  
  my @nvar = ();
  my @merr = ();
  my @verr = ();
  for my $n (2 .. $maxn)
  {
    print "$n";
    my ($m, $v, $merr, $verr) = Varsim($n, 30);
    push @nvar, $n;
    push @merr, $merr;
    push @verr, $verr;
    print "\n";
  }  
  
  my $nv = pdl(@nvar);
  my $em = pdl(@merr);
  my $ev = pdl(@verr);
  
  # linear regression of y = f / sqrt(n)
  my $fm = sum($em / sqrt($nv)) / sum(1 / $nv);
  my $fv = sum($ev / sqrt($nv-1)) / sum(1 / ($nv-1));
  print "F_mean = $fm\n";
  print "F_var = $fv\n";
  
  PlotPanelN($win, 1,
    BOX => [1, $maxn+1, 0, max($ev) * 1.1],
    XLAB => 'Number of replicates',
    YLAB => 'Fractional error',
    xbox => 'I', ybox => 'I',
    CHARSIZE => $cs
  );

  #LinePlot($win, $nv, $em, COLOR=>'BLUE');
  PointPlot($win, $nv, $em, COLOR=>'BLUE', SYMBOLSIZE=>1.5);
  #LinePlot($win, $nv, $ev, COLOR=>'RED');
  PointPlot($win, $nv, $ev, COLOR=>'RED', SYMBOLSIZE=>1.5);
  
  my $xx = sequence(1000)/10 + 2;
  my $ym = $fm / sqrt($xx);
  my $yv = $fv / sqrt($xx-1);
  LinePlot($win, $xx, $ym, COLOR=>'BLUE');
  LinePlot($win, $xx, $yv, COLOR=>'RED');
  
  FullBox($win);
  TopText($win, "Mean and variance median fractional error", c=>0.5, x=>0.5, y=>-2);
  TopText($win, $name, c=>0.5, x=>0.5, y=>-4);
  
  $win->close();
}

sub Varsim
{
  my ($nvar, $nsim) = @_;

  $nsim ||= 1000;
  
  my @var = ();
  my @mean = ();
  for my $k (1 .. $nsim)
  {
    #my $r = random($nrep);
    #my $is = qsorti $r;
    #my $rsel = $is(0:$nvar-1);
    
    my $sel = floor(random($nvar) * $nrep);  # sampling with replacement
    my ($ms, $ss) = statsover($d($sel,));
    push @mean, $ms;
    push @var, $ss**2;
  }

  my $mean = pdl(@mean)->transpose();
  my $var = pdl(@var)->transpose();
  my ($meanvar, $sdvar) = statsover($var);
  my ($meanmean, $sdmean) = statsover($mean);

  my $x = $lm;
  my $verr = $sdvar / $meanvar;
  my $merr = $sdmean / $meanmean;

  my ($dum1, $dum2, $m_verr) = stats($verr);
  my ($dum3, $dum4, $m_merr) = stats($merr);

  return $merr, $verr, $m_merr, $m_verr  
}

##############################################

sub PlotMeanVM
{
  my $win = NewWindow($psfile);

  my $vm = ($m / $s);  
  my ($miny, $maxy) = BoxMinMax($vm);
  #($miny, $maxy) = (-0.5, 4.5);
  
  
  PlotPanelN($win, 1,
    BOX => [$minm, $maxm, $miny, $maxy],
    XLAB => 'Mean',
    YLAB => 'M / S',
    xbox => 'L', ybox => ''
  );
  
  
  PointPlot($win, $lm, $vm, SYMBOLSIZE=>0.3);
  
  LinePlot($win, pdl($minm,$maxm), pdl(0,0), COLOR=>'RED', LINEWIDTH=>2);  # poisson
  Title($win);
  
  $win->close();
}

##############################################

sub PlotMeanlSM
{
  my $win = NewWindow($psfile);
  
  PlotPanelN($win, 1,
    BOX => [$minm, $maxm, 0, 1.5],
    XLAB => 'Mean (M)',
    YLAB => 'log(S) / log(M)',
    xbox => 'L', ybox => 'F'
  );
  
  my $x = $lm;
  my $y = $ls / $lm;

  #my $P = PoissonFromData($p);
  #my $sel = which($P < log10(0.05));
  #my $sel2 = which($P < log10(0.05/nelem($P)));

  
  PointPlot($win, $x, $y, SYMBOLSIZE=>0.3);
  #PointPlot($win, $x->dice($sel), $y->dice($sel), SYMBOLSIZE=>0.3);
  #PointPlot($win, $x->dice($sel2), $y->dice($sel2), SYMBOLSIZE=>0.5, COLOR=>'RED');
  
  my $xx = pdl(-2, 7);
  my $yy = pdl(0.5, 0.5);
  my $yy2 = pdl(1, 1);
  LinePlot($win, $xx, $yy, COLOR=>'RED', LINEWIDTH=>2);
  LinePlot($win, $xx, $yy2, COLOR=>'GOLD2', LINEWIDTH=>1);
  $win->text('Poisson', CHARSIZE=>0.8, COLOR=>'RED', TEXTPOSITION=>[5, 0.53, 0, 0, 0]);
  Title($win);
  
  $win->close();
}

##############################################

sub PlotMeanOver
{
  my ($log) = @_;
  
  my $win = NewWindow($psfile);
  
  my $x = log10($m);
  my $y = ($s*$s - $m) / ($m*$m);   # overdispersion
  my $mo = median($y);
  print "Median phi = $mo\n"; 
  
  my $sy = qsort $y;
  my $n = nelem($sy);
  my $ylo = at($sy, 0.05 * $n);
  my $yhi = at($sy, 0.95 * $n);
  
  if($log)
  {
    $y->where($y<=0) .= 1e-6;
    $y = log10($y);
  }
  
  my ($x1, $x2) = ($minm, $maxm);
  my ($y1, $y2) = ($log) ? (-6.2, 1) : (-0.3, 1.2);
  #($y1, $y2) = BoxMinMax($y);
  my %ybox = ($log) ? (ybox=>'L') : ();
  
  #($y1, $y2) = (-5, 1);
  
  PlotPanelN($win, 1,
    BOX => [$x1, $x2, $y1, $y2],
    XLAB => 'log mean',
    YLAB => '#gF',
    xbox => 'F', %ybox
  );
  
  PointPlot($win, $x, $y, SYMBOLSIZE=>0.3);
  LinePlot($win, pdl($x1,$x2), pdl(0,0), COLOR=>'RED') unless $log;
  
  LinePlot($win, pdl($x1,$x2), pdl($ylo,$ylo), COLOR=>'GOLD2') unless $log;
  LinePlot($win, pdl($x1,$x2), pdl($yhi,$yhi), COLOR=>'GOLD2') unless $log;
  
  PlotRunningMean($win, $x, $y, 100, 1);
  
  Title($win);
    
  $win->close();
  
}

##############################################

sub PlotDistOver
{
  my $win = NewWindow($psfile);
 
  my $o = ($s*$s - $m) / ($m*$m);   # overdispersion
  my $mean = mean($o);
  my $med = median($o);
  print "Median phi = $med\n";
  
  my ($min, $max) = (-0.02, 0.2);
  my $span = $max - $min;
  my $hmin = $min - 0.05 * $span;
  my $hmax = $max + 0.05 * $span;
  
  my ($x, $h) = BuildHistogram($o, 100, $hmin, $hmax, 1);
  
  #my ($x1, $x2) = BoxMinMax($x);
  my ($y1, $y2) = (0, max($h)*1.05);
  
  PlotPanelN($win, 1,
    BOX => [$min, $max, $y1, $y2],
    XLAB => 'Dispersion parameter',
    YLAB => 'Normalized count',
  );

  BarPlot($win, $x, $h, colour=>[215,215,255]);  
  
  Title($win);
    
  $win->close();
  
}


##############################################

sub CompareCumulVar
{
  my ($set1, $set2) = @_;
  
  my $win = NewWindow($psfile);
  
  PlotPanelN($win, 1,
    BOX => [0, 1, 0, 1],
    XLAB => 'S / M',
    YLAB => 'Fraction of genes',
    xbox => 'F', ybox => 'F'
  );
  
  my ($s1, $y1) = GetVar($set1);
  OverPlot($win, $s1, $y1, PLOTTYPE=>'LINE', COLOR=>'RED', LINEWIDTH=>2);
  my ($s2, $y2) = GetVar($set2);
  OverPlot($win, $s2, $y2, PLOTTYPE=>'LINE', COLOR=>'BLUE', LINEWIDTH=>2);
  
  $win->close();
}

##############################################

sub GetVar
{
  my $set = shift;
  
  my ($p, $m, $s) = ReadExpressionFile("");
  my $lm = log($m) / log(10);
  my $ls = log($s) / log(10);
  
  my $v = $ls / $lm;
  my $sv = qsort $v;
  my $N = nelem($sv);
  my $y = sequence($N) / $N;

  return ($sv, $y);  
}

##############################################

sub Composition
{
  my $win = NewWindow($psfile);
  
  $PL_xmin = 0.1;
  $PL_ymin = 0.1;
  PlotPanelN($win, 1,
    BOX => [0, 1, 0, 1],
    XLAB => 'Fraction of genes',
    YLAB => 'Fraction of total count',
    CHARSIZE => 0.9
  );

  my %cols = (WT => 'BLUE', Snf2 => 'RED');

  my $i = 0;
  for my $cond (@conds)
  {
    ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel) = @{$dat{$cond}};
    
    my $sm = qsort $m;
    my $n = nelem($sm);
    my @f = ();
    for my $i (0 .. $n-1)
    {
      my $x = at($sm, $n - 1 - $i);
      if($i == 0)
        {push @f, $x}
      else
         {push @f, $f[$i - 1] + $x}
    }
    
    LinePlot($win, sequence($n)/($n-1), pdl(@f)/sum($sm), COLOR=>$cols{$cond}, LINEWIDTH=>3);

    my $x0 = 0.7;
    my $y0 = 0.4;
    my $dx = 0.05;
    my $dy = 0.05;
    my $y = $y0 - $i * $dy;
    LinePlot($win, pdl($x0,$x0+$dx), pdl($y,$y), LINEWIDTH=>3, COLOR=>$cols{$cond});
    $win->text($cond, COLOR=>$cols{$cond},  TEXTPOSITION=>[$x0+$dx+0.02,$y, 0, 0, 0]);
    
    $i++
  }
  
  
  $win->close();
}
    
##############################################

sub sellist
{
  my ($g2i) = @_;
  my @hl = ();
  if(defined $sellist)
  {
    my @selfiles = split(/\,/, $sellist);
    for my $selfile (@selfiles)
    {
      die "Cannot find file $selfile\n" unless -e $selfile;
      my @selgen = ReadGeneList($selfile);
      my @s = map {$g2i->{$_}} @selgen;
      push @hl, pdl(@s);
    }
  }
  return @hl;
}



sub Title
{
  my ($win) = @_;
  FullBox($win);
  TopText($win, $name, x=>0.07, y=>-2, CHARSIZE => 0.7);
  TopText($win, "Normalization: $norm", x=>0.07, y=>-4, CHARSIZE => 0.7);
}

##############################################

sub PlotLoess
{
  my ($w, $x, $y) = @_;
  
  return unless $withloess;
  
  print "Calculating loess...";
  my ($xs, $ys) = Loess($x, $y);
  print " done\n";
  
  LinePlot($w, $xs, $ys, COLOR=>[110,110,255], LINEWIDTH=>2);
  
  
}

##############################################

sub Loess
#
# Smooth profiles using R function loess
#
{
  my ($x, $y, $span, $degree) = @_;
  $span = 0.75 if !defined $span;
  $degree = 2 if !defined $degree;
  
  my $i = qsorti $x;
  my $xs = $x($i);
  my $ys = $y($i);
  
  my $xyfile = "/tmp/xy$$.dat"; 
  my $smoothfile = "/tmp/loess$$.dat";
  my $Rfile = "/tmp/loess$$.R";
  unlink $xyfile, $smoothfile, $Rfile;
  wcols $xs, $ys, $xyfile;

  local *R;
  open R, ">$Rfile" or die "Cannot open temporary file $Rfile\n";
  print R <<EOF;
dat = read.table("$xyfile")
m = as.matrix(dat)
x = m[,1]
y = m[,2]
ls = loess(y ~ x, span=$span, degree=$degree)
pr = predict(ls, x)
write.table(pr, file="$smoothfile", row.names=FALSE, col.names=FALSE)
EOF
  close R;
  call("R --vanilla < $Rfile", my $stat);
  die "Problem with R call Loess\n"if !-e $smoothfile || $stat;

  my $smooth = rcols $smoothfile; 

  unlink $xyfile, $smoothfile, $Rfile;
  return $xs, $smooth;
}


=head1 SYNOPSIS

  plot_ms.pl
  plot_ms.pl -clean -withloess -psfile=figure.ps

=head1 OPTIONS

=over 4

=item B<-wthloess>

If specified, loess fit will be added to the plot.

=item B<-norm>=I<string>

Expression data normalization. Input expression files contain raw counts. These are normalized on the fly, using one of the following methods:

  none - no normalization
  lev - levelling (equal read count) normalization (needs data preparation)
  deseq - DESeq normalization (default)
  tmm - trimmed mean of M values normalization
  fpkm - approximate FPKM normalization (not identical to cuffdiff!)
  totcount - total count
  totcountm - total count with normalization factors scaled to the mean of 1

=item B<-clean>

Use clean data replicates, as specified in C<defs.dat> file.

=item B<-psfile>=I<name>

Name of the postsript file to redirect output to.

=back
