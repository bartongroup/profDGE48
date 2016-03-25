#!/sw/bin/perl

=head1 NAME

B<compare_powerstats_db.pl>

=head1 DESCRIPTION

Script to create figure 2.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use Math::CDF qw(:all);

use GRNASeq;
use GetOpt;
use Powertests;

use CompBio::PLgraphs;
use CompBio::Tools;
use CompBio::HTMLTable;
use CompBio::Statistics;
use CompBio::Distribution;
use CompBio::Statistics;

use Switch 'Perl6';


$| = 1;
Stamp();

#my $multicor = 'bh';

my $testinfofile = "de_tests.txt";
my $dir = 'powerstats_db_ref';
my $truetest = 'self';
my $tests = 'lt,bayseq,cuffdiff,degseq,deseq1,deseq2,ebseq,edger,limma,noiseq,samseq,poissonseq';
my $minnrep = 3;
my $witherrors;
my $what = 'tpr';
my $fclim = 0.5;
my $nrep = 3;
my $graph = 'rates';
my ($help, $man);
GetMyOptions(
  'dir=s' => \$dir,
  'testinfofile=s' => \$testinfofile,
  'truetest=s' => \$truetest,
  'tests=s' => \$tests,
  'minnrep=i' => \$minnrep,
  #'what=s' => \$what,  # do not change
  #'fclim=f' => \$fclim,  # not in use
  'nrep=i' => \$nrep,
  'graph=s' => \$graph,
  #witherrors => \$witherrors,  # not in use
   help => \$help,
   man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;

my $cs = 0.45;
my $cs2 = 0.3;
my @fclim = (0.1, 0.3, 0.5, 1, 2);
my @nrep = (3, 5, 10, 20, 30);
my @col = @contrast_colours;

@nrep = (2, @nrep) if $minnrep == 2;

###############################################

my $maxnrep = 42;


my %tests = ReadTestfileInfo($testinfofile);
my @tests = reverse split(/\,/, $tests);
my $N = @tests;

# Plot power curves

my $div = 0.8;



$PL_nx = 3;
$PL_ny = 3;
$PL_deltax = 0.07;
$PL_deltay = 0.07;
$PL_xmax = $div;

my $w = NewWindow($psfile);



given($graph)
{
  when 'rates' {CompareRatesCombined()}
  when 'rocs' {PlotROCs()}
  when 'sp' {PlotSigProps()}
  when 'ns' {PlotNsigs()}
  default {die "Unknown -graph\n"}
}

Legend($w) unless $graph eq 'rates';


#FullBox($w);
#TopText($w, "", c=>0.5, x=>0.35, y=>-2);

$w->close();


sub PlotNsigs
{
  my $pan = 1;
  
  $PL_deltax = 0.03;
  $PL_nx = 3;
  
  PlotNsig($w, $pan++, 0);
  PlotNsig($w, $pan++, 0.5);
  #PlotNsig($w, $pan++, 1);
  PlotNsig($w, $pan++, 2);
}


sub PlotROCs
{
  my $pan = 1;
  
  PlotROC($w, $pan++, 3, 0);
  PlotROC($w, $pan++, 3, 0.5);
  PlotROC($w, $pan++, 3, 1);
  
  PlotROC($w, $pan++, 6, 0);
  PlotROC($w, $pan++, 6, 0.5);
  PlotROC($w, $pan++, 6, 1);

  PlotROC($w, $pan++, 12, 0);
  PlotROC($w, $pan++, 12, 0.5);
  PlotROC($w, $pan++, 12, 1);
}

sub PlotSigProps
{
  my $pan = 1;
  
  $PL_deltax = 0.04;
  
 PlotSigProp($w, $pan++, 3); 
 PlotSigProp($w, $pan++, 6);
 PlotSigProp($w, $pan++, 20);
}

sub PlotRates
{
  my $pan = 1;
  #my $what = 'fpr';
  
  PlotNrepRateFC($w, $pan++, 0, $what);
  PlotNrepRateFC($w, $pan++, 0.1, $what);
  PlotNrepRateFC($w, $pan++, 0.3, $what);
  PlotNrepRateFC($w, $pan++, 0.5, $what);
  PlotNrepRateFC($w, $pan++, 1, $what);
  PlotNrepRateFC($w, $pan++, 2, $what);
}

sub CompareRates1
{
  my $pan = 1;
  #my $what = 'p';
  #my $fclim = 0.5;
  
  CompareRateFC($w, $pan++, 3, $fclim, $what);
  CompareRateFC($w, $pan++, 4, $fclim, $what);
  CompareRateFC($w, $pan++, 5, $fclim, $what);
  CompareRateFC($w, $pan++, 6, $fclim, $what);
  CompareRateFC($w, $pan++, 8, $fclim, $what);
  CompareRateFC($w, $pan++, 12, $fclim, $what);
  CompareRateFC($w, $pan++, 16, $fclim, $what);
  CompareRateFC($w, $pan++, 20, $fclim, $what);
  CompareRateFC($w, $pan++, 30, $fclim, $what);
}
    
sub CompareRates2
{
  my $pan = 1;
  #my $what = 'fpr';
  #my $nrep = 3;
  
  CompareRateFC($w, $pan++, $nrep, 0, $what);
  #CompareRateFC($w, $pan++, $nrep, 0.3, $what);
  CompareRateFC($w, $pan++, $nrep, 0.5, $what);
  #CompareRateFC($w, $pan++, $nrep, 0.8, $what);
  #CompareRateFC($w, $pan++, $nrep, 1, $what);
  CompareRateFC($w, $pan++, $nrep, 2, $what);
}

sub CompareRates3
{
  my $pan = 1;
  #my $what = 'fpr';
  #my $nrep = 3;
  
  $what = 'tpr';
  CompareRateFC($w, $pan++, $nrep, 0, $what);
  CompareRateFC($w, $pan++, $nrep, 0.5, $what);
  CompareRateFC($w, $pan++, $nrep, 2, $what);
  $what = 'fpr';
  CompareRateFC($w, $pan++, $nrep, 0, $what);
  CompareRateFC($w, $pan++, $nrep, 0.5, $what);
  CompareRateFC($w, $pan++, $nrep, 2, $what);
  
  FullBox($w);
  TopText($w, "Nrep = $nrep", c=>0.5, x=>0.35, y=>-2);
}

sub CompareRatesCombined
{
  $PL_ymin = 0.65;
  $PL_ymax = 0.95;
  $div = 0.3;
  Names();
  
  $PL_xmin = $div;
  $PL_xmax = 0.95;
  $PL_nx = 3;
  $PL_deltax = 0.03;
  
  my $pan = 1;
  CompareRateTrueFalse($w, $pan++, $nrep, 0); 
  CompareRateTrueFalse($w, $pan++, $nrep, 0.5);
  CompareRateTrueFalse($w, $pan++, $nrep, 2);
}

##########################################

sub PanelNrep01
{
  my ($w, $pan, $lab, $min, $max) = @_;
  
  $min ||= 0;
  $max ||= 1;
  PlotPanelN($w, $pan,
    BOX => [0, $maxnrep, $min, $max],
    XLAB => 'Number of replicates',
    YLAB => $lab,
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    #forceylab => 1,
    YTICK => 0, NYSUB => 0,
    CHARSIZE => $cs
  );
}

sub PanelFC01
{
  my ($w, $pan, $lab) = @_;
  
  PlotPanelN($w, $pan,
    BOX => [0, 3, 0, 1],
    XLAB => 'log#d2#uFC',
    YLAB => $lab,
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    forceylab => 1,
    YTICK => 0, NYSUB => 0,
    CHARSIZE => $cs
  );
}
  
##########################################

sub CompareRateFC
{
  my ($w, $pan, $n, $fc, $what) = @_;
  
  $what ||= 'fpr';
  my %lab = (tpr => 'True positive rate', fpr => 'False positive rate');
  
  my $max = 1;
  $max = 0.3 if $what eq 'fpr';
  $max = 0.07 if $what eq 'fdr';
  $max = 6000 if $what eq 'p';
  PlotPanelN($w, $pan,
    BOX => [0, $N+1, 0, $max],
    XLAB => '',
    YLAB => $lab{$what},
    XBOX =>'B',
    ybox => 'NI',
    CHARSIZE => $cs
  );
  
  my ($lab, $r, $sr);
  for my $i (0 .. $N - 1)
  {
    my $test = $tests[$i];
    my $tt = ($truetest eq 'self') ? $test : $truetest;
    my $fcfile2 = "$dir/power_${test}_true-${tt}_fc2.stats";
    
    my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetFCData($fcfile2);
    if($what eq 'fnr')
    {
      $lab = 'FNR';
      ($r, $sr) = rerr($fn, $tp, $sfn, $stp)
    }
    elsif($what eq 'fpr')
    {
      $lab = 'FPR';
      ($r, $sr) = rerr($fp, $tn, $sfp, $stn);
    }
    elsif($what eq 'tpr')
    {
      $lab = 'TPR';
      ($r, $sr) = rerr($tp, $fn, $stp, $sfn)
    }
    elsif($what eq 'p')
    {
      $lab = 'P';
      ($r, $sr) = ($tp+$fp, sqrt($stp**2 + $sfp**2))
    }
    elsif($what eq 'fdr')
    {
      $lab = 'FDR';
      ($r, $sr) = rerr($fp, $tp, $sfp, $stp)
    }
    
    #my ($r, $sr) = rerr($tp, $fn, $stp, $sfn);
    my $sel = which(($fclim == $fc) & ($nrep == $n));
    my ($rf, $srf) = ($r($sel), $sr($sel));
    #print "$pan:  $rf $srf  ",  $tp($sel), ":", $tn($sel), ":", $fp($sel), ":", $fn($sel), "\n";
    
    PointPlot($w, pdl($i+1), $rf, YERRORBAR=>2*$srf, SYMBOLSIZE=>1, COLOR=>$col[$i]);
  }
  #TopText($w, "|lg#d2#uFC|>$fc N#drep#u=$n", y=>1, CHARSIZE=>$cs);
  TopText($w, "|lg#d2#uFC|>$fc", y=>1, CHARSIZE=>$cs);
}


##########################################

sub CompareRateTrueFalse
{
  my ($w, $pan, $n, $fc) = @_;
  
  #my %lab = (tpr => 'True positive rate', fpr => 'False positive rate');
  
  tfpanel($w, $pan, $n, $fc, 'fpr');
  tfpanel($w, $pan, $n, $fc, 'tpr');

  PlotPanelN($w, $pan, %PL_empty);
  TopText($w, "|lg#d2#uFC|>$fc", y=>2, x=>0.5, c=>0.5,CHARSIZE=>$cs);
}

sub tfpanel
{
  my ($w, $pan, $n, $fc, $what) = @_;
  
  my %lab = (tpr => 'TPR', fpr => 'FPR');
  my $BAD = -666;
   
  # panel coordinates
  my $x = ($pan - 1) % $PL_nx + 1;
  my $y = int(($pan - 1) / $PL_nx) + 1;
  my ($x1, $x2, $y1, $y2) = PanelCoord($x, $y);
  my $xdiv = ($x2 + $x1) / 2;

  my ($p1, $p2) = ($what eq 'tpr') ? ($xdiv, $x2) : ($x1, $xdiv);
  #my ($xtick, $nxsub) = ($what eq 'tpr') ? (0.5, 5) : (0.2, 2);
  my ($xtick, $nxsub) = ($what eq 'tpr') ? (0.5, 5) : (0.5, 5);

  #my ($min, $max) = ($what eq 'tpr') ? (0, 1) : (-0.6, 0);
  my ($min, $max) = ($what eq 'tpr') ? (0, 1) : (-1, 0);
  $w->xyplot(pdl(-1), pdl(-1),
    VIEWPORT => [$p1, $p2, $y1, $y2],
    BOX => [$min, $max, 0, $N+1],
    %PL_empty,
    XBOX=>'BC', YBOX=>'BC'
  );
  
  
  my @x = ();
  my @y = ();
  my @sy = ();
  for my $i (0 .. $N - 1)
  {
    my $test = $tests[$i];
    my $tt = ($truetest eq 'self') ? $test : $truetest;
    my $fcfile2 = "$dir/power_${test}_true-${tt}_fc2.stats";
    my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetFCData($fcfile2);
    my ($r, $sr);
    if($what eq 'fpr')
      {($r, $sr) = rerr($fp, $tn, $sfp, $stn)}
    elsif($what eq 'tpr')
      {($r, $sr) = rerr($tp, $fn, $stp, $sfn)}
    
    if($what eq 'fpr')  # left panel
    {
      $r *= -1;
      $sr *= -1
    }
    
    my $sel = which(($fclim == $fc) & ($nrep == $n));
    my ($rf, $srf) = ($r($sel), $sr($sel));
    push @x, $i + 1;
    my $_rf = at($rf, 0);
    my $_srf = at($srf, 0);
    $_rf = $BAD if $_rf =~/nan/;
    $_srf = $BAD if $_srf =~/nan/;
    push @y, $_rf;
    push @sy, $_srf;
  }
  my $xx = pdl(@x);
  my $yy = pdl(@y);
  my $syy = pdl(@sy);
  $yy->badflag(1);$yy->badvalue($BAD);
  $syy->badflag(1);$syy->badvalue($BAD);
  my %o = (y0 => 1) if $what eq 'fpr';
  my $c = ($what eq 'fpr') ? [215, 215, 255] : [215, 215, 165];
  my $cd = ($what eq 'fpr') ? [15, 15, 215] : [110, 110, 25];
  #my @cd = map {$_ * 0.5} @$c;
  #print join(":", @cd), "\n";
  
  BarPlot($w, $xx, $yy, colour=>$c, swap=>1, boxes=>1, width => 1, boxlinewidth=>1, boxlinecolour=>[120,120,120], %o);
  PointPlot($w, $yy, $xx, XERRORBAR=>2*$syy, SYMBOLSIZE=>0.9, COLOR=>$cd);

   vgridlines($w, $min, $max, $xtick);
  $min *= -1;
  my $ybox = ($what eq 'tpr') ? 'B' : '';
  $w->xyplot(pdl(-1), pdl(-1),
    VIEWPORT => [$p1, $p2, $y1, $y2],
    BOX => [$min, $max, 0, $N+1],
    YLAB => '',
    XLAB => $lab{$what},
    YBOX =>$ybox,
    XBOX => 'BNSTI',
    XTICK => $xtick,
    NXSUB => $nxsub,
    CHARSIZE => $cs
  );
  
}

sub vgridlines
{
  my ($w, $x1, $x2, $xs) = @_;
  
  my $col = [120,120,120];
  for(my $x = $x1; $x <= $x2; $x += $xs)
    {LinePlot($w, pdl($x,$x), pdl(-1e4,1e4), COLOR=>$col)}
}


##########################################

sub PlotNrepRateFC
{
  my ($w, $pan, $fc, $what) = @_;
  $what ||= 'tpr';
  
  my $min = 0.8 if $what eq 'tnr';
  my $max = 0.2 if $what eq 'fpr';
  PanelNrep01($w, $pan, $what, $min, $max);
  my $i = 0;
  for my $test (@tests)
  {
    LineNrepRateFC($w, $test, $fc, $col[$i], $what);
    $i++
  }
  TopText($w, "|lg#d2#uFC|>$fc", x=>0.95, c=>1, y=>2,  CHARSIZE=>0.4);
}  

##########################################

sub PlotSigProp
{
  my ($w, $pan, $nrep) = @_;
  
  print "$nrep:  ";
  PlotPanelN($w, $pan,
    BOX => [0, 1, 0, 1], %PL_empty
  );
  
  TopText($w, $nrep, y=>1);
  
  for my $i (0 .. @tests - 1)
  {
    my $test = $tests[$i];
    print "$test ";
    my $file = "./powerstats_db_ref/power_${test}_true-${test}_sigprop.stats";
    next unless -e $file;
    
    my ($n, $prop) = rcols $file, 0, 3;
    my $sel = which($n == $nrep);
    my $p = $prop($sel);
    my $ps = qsort $p;
    my $x = (sequence($ps) + 1) / nelem($ps);
  
    LinePlot($w, $x, $ps, LINEWIDTH=>3, COLOR=>$contrast_colours[$i]);
  }
  print "\n";
  
  PlotPanelN($w, $pan,
    BOX => [0, 1, 0, 1],
    YLAB => 'SDE bootstrap proportion',
    XLAB => 'Gene fraction',
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    #forceylab => 1,
    CHARSIZE => $cs
  );
  PlotPanelN($w, $pan,
    BOX => [0, 1, 0, 1],
    XBOX => 'C', YBOX => 'C',
  );
  
}

##########################################

sub PlotNsig
{
  my ($w, $pan, $fc) = @_;

  my ($min, $max) = (0, 7000);
  $max = 2000 if $fc == 0.5;
  $max = 700 if $fc == 1;
  $max = 150 if $fc == 2;
  
  PanelNrep01($w, $pan, 'Numbe of SDE genes', $min, $max);
  my $i = 0;
  for my $test (@tests)
  {
    LineNrepRateFC($w, $test, $fc, $col[$i], 'nsig');
    $i++
  }
  TopText($w, "|lg#d2#uFC|>$fc", x=>0.5, c=>0.5, y=>2,  CHARSIZE=>0.4);
}  

##########################################

sub CompareNsigNull
{
  #my @tests = split(/\,/, $tests);
  #my $N = @tests;
  #my %tests = ReadTestfileInfo($testinfofile);
  
  my $div = 0.3;
  
  $PL_xmin = 0.05;
  $PL_xmax = $div;
  PlotPanelN($w, 1, %PL_empty, BOX => [0,1,0,$N+1]);
  for my $i (1 .. $N)
  {
    my $name = $tests{$tests[$i-1]}{name};
    $w->text($name, COLOR=>'BLACK', CHARSIZE=>0.7, TEXTPOSITION=>[0.9, $i, 0, 0, 1])
  }
    
    $PL_xmin = $div;
    $PL_xmax = 0.95;
    PlotPanelN($w, 1,
      BOX => [0, 0.2, -1, $N],
      XBOX => 'BNST',
      YBOX => 'B',
      XLAB => 'Fraction of significant genes'
    );
    
  my @dat = ();
  local *F;
  
  for my $test (@tests)
  {
    my $file = "powerstats_same_db/power_${test}_true-${test}_nsig.stats";
    open F, $file or die "Cannot open $file\n";
    <F>;
    while(my $line = <F>)
    {
      my @d = split /\t/, $line;
      my $n = shift @d;
      next unless $n == $nrep;
      my $d = pdl(@d);
      push @dat, $d;
    }
  }
    
  BoxPlot($w, \@dat, cloud=>1, swap=>1, pcolour=>[110,110,110]);
  LinePlot($w, pdl(0.05,0.05), pdl(-1,$N), COLOR=>'RED');
}

sub errorbar
{
  my ($win, $xx, $yy, $ee1, $ee2, $col) = @_;

  $col ||= 'BLACK';
  my $w = 1;

  for my $i (0 .. nelem($xx) - 1)
  {
    my $x = at($xx, $i);
    my $y = at($yy, $i);
    my $e1 = at($ee1, $i);
    my $e2 = at($ee2, $i);
    my $s = 0.03;
    my $xp = pdl($x-$s, $x+$s, $x, $x, $x-$s, $x+$s);
    my $yp = pdl($y-$e1, $y-$e1, $y-$e1, $y+$e2, $y+$e2, $y+$e2);
    OverPlot($win, $xp, $yp, PLOTTYPE=>'LINE', LINEWIDTH=>$w, COLOR=>$col, LINESTYLE=>1);
  }
}

##########################################

sub LineNrepRateFC
{ 
  my ($w, $test, $fc, $col, $what) = @_;
  
  $what ||= 'tpr';
  my $tt = ($truetest eq 'self') ? $test : $truetest;
  my $fcfile2 = "$dir/power_${test}_true-${tt}_fc2.stats";
  
  my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetFCData($fcfile2);
  selrep($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn);
  my $maxn = max($nrep);

  my ($lab, $r, $sr);
  if($what eq 'fnr')
  {
    $lab = 'FNR';
    ($r, $sr) = rerr($fn, $tp, $sfn, $stp)
  }
  elsif($what eq 'fpr')
  {
    $lab = 'FPR';
    ($r, $sr) = rerr($fp, $tn, $sfp, $stn)
  }
  elsif($what eq 'tpr')
  {
    $lab = 'TPR';
    ($r, $sr) = rerr($tp, $fn, $stp, $sfn)
  }
  elsif($what eq 'tnr')
  {
    $lab = 'TNR';
    ($r, $sr) = rerr($tn, $fp, $stn, $sfp)
  }
  elsif($what eq 'nsig')
  {
    $lab = 'N#dsig#u';
    ($r, $sr) = ($tp + $fp, sqrt($stp**2 + $stn**2));
  }


#wcols $nrep, $fclim, $fn, $tp, $r, $sr;die;
  
  my @n = (); my @R = (); my @sR = ();
  for my $n (2 .. $maxn)
  {
    my $sel = which(($fclim == $fc) & ($nrep == $n));
    next if nelem($sel) == 0;
      #print "S = $sel\n";
    my ($rf, $srf) = ($r($sel), $sr($sel));
    push @n, $n;
    push @R, sum($rf);
    push @sR, sqrt(sum($srf**2));
  }
  if($witherrors)
    {linewitherr($w, pdl(@n), pdl(@R), pdl(@sR), $col)}
  else
    {LinePlot($w, pdl(@n), pdl(@R), COLOR=>$col, LINEWIDTH=>3)}
    #wcols pdl(@n), pdl(@R), pdl(@sR);die;
}

##########################################

sub LineFCRateNrep
{ 
  my ($w, $test, $n, $col, $what) = @_;
  
  $what ||= 'TPR';
  my $tt = ($truetest eq 'self') ? $test : $truetest;
  my $fcfile2 = "$dir/power_${test}_true-${tt}_fc2.stats";

  my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetFCData($fcfile2);
  selrep($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn);
  #my $maxn = max($nrep);
  # collect all fclim values
  my @F = CollectFC($fclim);

  my ($lab, $r, $sr);
  if($what eq 'fnr')
  {
    $lab = 'FNR';
    ($r, $sr) = rerr($fn, $tp, $sfn, $stp)
  }
  elsif($what eq 'fpr')
  {
    $lab = 'FPR';
    ($r, $sr) = rerr($fp, $tn, $sfp, $stn)
  }
  elsif($what eq 'tpr')
  {
    $lab = 'TPR';
    ($r, $sr) = rerr($tp, $fn, $stp, $sfn)
  }
  

  #wcols $nrep, $fclim, $fn, $tp, $r, $sr;die;

  
  my @f = (); my @R = (); my @sR = ();
  for my $f (@F)
  {
    my $sel = which(($fclim == $f) & ($nrep == $n));
    next if nelem($sel) == 0;
      #print "S = $sel\n";
    my ($rf, $srf) = ($r($sel), $sr($sel));
    push @f, $f;
    push @R, sum($rf);
    push @sR, sqrt(sum($srf**2));
  }
  linewitherr($w, pdl(@f), pdl(@R), pdl(@sR), $col);
}

#######################################


sub PlotROC
{
  my ($w, $pan, $nroc, $froc) = @_;
  
  PlotPanelN($w, $pan,
    BOX => [0, 0.2, 0, 1],
    XLAB => 'FPR',
    YLAB => 'TPR',
    xbox =>'NI',
    ybox => 'NI',
    CHARSIZE => $cs
  );
  LinePlot($w, pdl(0,1), pdl(0,1), COLOR=>'RED', LINESTYLE=>4);
  LinePlot($w, pdl(-1,-1), pdl(-1,-1), LINESTYLE=>1);
  
  my @fpr_err = ();
  my @tpr_err = ();
  for my $i (0 .. $N - 1)
  {
    my $test = $tests[$i];
    my $tt = ($truetest eq 'self') ? $test : $truetest;
    my $rocfile = "$dir/power_${test}_true-${tt}_roc.stats";

    my ($nrep, $fclim, $fpr, $tpr, $sfpr, $stpr) = rcols $rocfile, {LINES=>'1:-1'};
    my $sel = which(($nrep == $nroc) & ($fclim == $froc) & ($fpr > 0) & ($tpr > 0));
    $fpr = $fpr($sel);
    $tpr = $tpr($sel);
    $sfpr = $sfpr($sel);
    $stpr = $stpr($sel);
    #print "$sfpr\n";
    LinePlot($w, $fpr, $tpr, COLOR=>$col[$i], LINEWIDTH=>3);
    my $esel = which(($fpr>0.1) & ($tpr>0.6));  # for typical errors only
    push @fpr_err, $sfpr($esel);
    push @tpr_err, $stpr($esel);
  }
  my ($fm, $fs, $fpr_err) = stats(pdl(@fpr_err)->flat()); 
  my ($tm, $ts, $tpr_err) = stats(pdl(@tpr_err)->flat());
  PointPlot($w, pdl(0.7), pdl(0.3), XERRORBAR=>2*$fpr_err, YERRORBAR=>2*$tpr_err, COLOR=>'GREY');
  TopText($w, "|lg#d2#uFC|>$froc N#drep#u=$nroc", y=>1, CHARSIZE=>$cs);
}

sub Legend
{
  $PL_xmin = $div;
  $PL_xmax = 0.99;
  $PL_ymax = 0.99;
  $PL_nx = 1;
  $PL_ny = 1;
  PlotPanelN($w, 1, %PL_empty);
  
  my $i = 0;
  for my $i (0 .. @tests - 1)
  {
    my $x = 0.01;
    my $y = 0.99 - $i * 0.04;
    my $dx = 0.18;
    
    my $c = $col[$i];
 
    LinePlot($w, pdl($x, $x+$dx), pdl($y, $y), COLOR=>$c, LINEWIDTH=>3);
    $w->text($tests{$tests[$i]}{name}, COLOR=>$c, CHARSIZE=>0.5, TEXTPOSITION=>[$x+$dx+0.05,  $y, 0, 0, 0]);
  }
}  


sub Names
{
  $PL_xmin = 0.05;
  $PL_xmax = $div;
  $PL_nx = 1;
  $PL_ny = 1;
  $PL_deltax = 0;
  
  my $N = @tests;
  
  PlotPanelN($w, 1, %PL_empty, BOX=>[0, 1, 0, $N+1]);
  
  my $i = 0;
  for my $i (0 .. @tests - 1)
  {
    my $x = 0.95;
    my $y = $i + 1;
    
    #my $c = $col[$i];
    my $c = 'BLACK';
    $w->text($tests{$tests[$i]}{name}, COLOR=>$c, CHARSIZE=>0.5, TEXTPOSITION=>[$x, $y, 0, 0, 1]);
  }
}  

  
sub linewitherr
{
  my ($w, $x, $r, $sr, $col) = @_;
  LinePlot($w, $x, $r, LINEWIDTH=>3, COLOR=>$col);
  LinePlot($w, $x, $r+$sr, LINEWIDTH=>1, COLOR=>$col);
  LinePlot($w, $x, $r-$sr, LINEWIDTH=>1, COLOR=>$col);
}

sub selrep
{
  my $nrep = $_[0];
  my $sel = which($nrep >= $minnrep);
  for my $i (0 .. @_ - 1)
    {$_[$i] = $_[$i]->($sel)}
}

=head1 SYNOPSIS

  compare_powerstats_db.pl     
    
=head1 OPTIONS

=item B<-testinfofile>=I<pathname>

A file with DE tools metadata. See C<plot_fcp.pl> for details. Default name is C<de_tests.txt>

=item B<-dir>=I<string>

Directory where the results from C<make_powerstats_db.pl> are stored. The default value is C<powerstats_db_ref>.


=item B<-truetest>=I<string>

A test to be used as a 'gold standard'. The default value (as used in the Paper) is 'self', that is the full clean replicate test for the same tool.

=item B<-tests>=I<string>

A comma-delimited list of tests/tools to be included in the plot. The default value is 'lt,bayseq,cuffdiff,degseq,deseq1,deseq2,ebseq,edger,limma,noiseq,samseq,poissonseq'.

=item B<-minnerp>=I<integer>

Minumum number of replicates to include in the plots. Default is 3.

=item B<-nrep>=I<number>

Number of replicates to create the plot.


=item B<-graph>=I<string>

Type of graph to plot. The default value is 'rates'. Do not change!


=cut
