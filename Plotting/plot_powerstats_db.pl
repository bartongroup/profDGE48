#!/sw/bin/perl

=head1 NAME

B<plot_powerstats_de.pl>

=head1 DESCRIPTION

Plot various results.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use Math::CDF qw(:all);
use Switch 'Perl6';

use GRNASeq;
use GetOpt;
use Powertests;

use PLgraphs;
use Tools;
use HTMLTable;
use Stats;
use Distribution;

my $version=1.0;

$| = 1;
Stamp();


#my $multicor = 'bh';

#my $tests = 't,lt,bayseq,cuffdiff,deseq,edger,limma,noiseq,samseq,poissonseq';

my $root;
my $title;
my $test;
my $nroc = 3;
my $froc= 0;
my $minnrep = 3;
my $maxfc = 2;
my $nrep = 3;
my $page = 1;
my $testinfofile = './de_tests.txt';
my $powerdir = 'powerstats_db_ref';
my $graph = 'over';
my $cs = 0.45;
my $cs2 = 0.3;
my ($help, $man);
GetMyOptions(
  #'root=s' => \$root, #obsolete
  'title=s' => \$title,
  'test=s' => \$test,
  #'nroc=i' => \$nroc, # not in use
  #'froc=f' => \$froc, # not in use
  #'page=i' => \$page, # not in use
  'testinfofile=s' => \$testinfofile,
  'powerdir=s' => \$powerdir,
  'minnrep=i' => \$minnrep,
  #'nrep=i' => \$nrep, # not in use
  'graph=s' => \$graph,
  'maxfc=f' => \$maxfc,
  'cs=f' => \$cs,
  'cs2=f' => \$cs2,
   help => \$help,
   man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;

die "Need -root or -test\n" unless defined $root || defined $test;

$root ||= "$powerdir/power_${test}_true-${test}";
print "Using root $root\n";

my $powerfile = "$root.stats";
my $fcfile1 = "${root}_fc1.stats";
my $fcfile2 = "${root}_fc2.stats";
my $rocfile = "${root}_roc.stats";
my $nsigfile = "${root}_nsig.stats";
my $sigpropfile = "${root}_sigprop.stats";

my @fclim = (0, 0.3, 1, 2);
my @nrep = (3, 6, 10, 20, 30);
#my @col = qw(BLACK BLUE RED GOLD2 GREEN MAGENTA GREY);
my @col = @contrast_colours;

@nrep = (2, @nrep) if $minnrep == 2;

###############################################


my $w;

given($graph)
{
  when 'over' {Overview($page)}
  when 'sim' {Simple()}
  when 'cn' {CompareNsig()}
  when 'ns' {PlotNsig()}
  when 'sp' {SigProp()}
  default {die "Unknown -graph\n"}
}


##########################################

sub Simple
{
  $PL_nx = 2;
  $PL_ny = 1;
  $PL_ymin = 0.5;
  $PL_deltax = 0.07;
  $cs = 0.7;
  
  $w = NewWindow($psfile);
  
  my $pan = 1;
  PlotNrepRateFC($pan++, 'tpr', maxn=>40);
  PlotFCRateNrep($pan++, 'tpr');
  
  if($root =~ /\/power_(\w+)_/)
  {
    my %tests = ReadTestfileInfo($testinfofile);
    my $name = $tests{$1}{name};
    if(defined $name)
    {
      FullBox($w);
      TopText($w, $name, CHARSIZE=>1.0, x=>0.5, c=>0.5, y=>-1);
    }
  }
  

  $w->close();  
}

##########################################

sub Overview
{
  my ($page) = @_;
  
  $page ||= 1;
  
  $PL_nx = 2;
  $PL_ny = 2;
  $PL_ymin = 0.05;
  $PL_deltax = 0.06;
  $PL_deltay = 0.07;
  $PL_xmin = 0.08;
  $PL_ymin = 0.08;
  $PL_xmax = 0.99;
  $cs = 0.5;
  $cs2 = 0.35;
  
  $w = NewWindow($psfile);
  
  my $pan = 1;
  
  PlotNsig($pan++);
  #PlotNrepRate2($pan++);
  PlotNrepRateFC($pan++, 'tpr', withfalse=>1);
  #PlotFCRateNrep($pan++, 'tpr', withfalse=>1);

  my $first = shift @nrep if $nrep[0] == 2;
  PlotFCRateNrep($pan++, 'tpr', withfalse=>1);
  @nrep = (2, @nrep) if $minnrep == 2;

  #PlotSigPropCumul($pan++, 16);
  PlotFour($pan++);
  #PlotFour($pan++, 1);
  #PlotROC($pan++);
  #PlotFC(3, 'fnr');
  #PlotFC($pan++, 'tpr');
  #PlotNrepRateFC(4, 'sigfrac');
  #PlotNrepRateFC($pan++, 'fpr', ylegend=>0.9);
  #PlotFCRateNrep(6, 'sigfrac');
  #PlotNrep($pan++, $nrep, $tp, $stp, 'TP');
  #PlotNrep(2, $nrep, $tn, $stn, 'TN');
  #PlotNrep(3, $nrep, $fp, $sfp, 'FP');
  #PlotNrep(4, $nrep, $fn, $sfn, 'FN');
  #PlotNrep(5, $nrep, $fpr, $sfpr, 'FPR');
  
  NameTitle($w);  
  $w->close();
}

##########################################

sub NameTitle
{
  my $w = shift;
  if(!defined $title && $root =~ /\/power_(\w+)_/)
  {
    my %tests = ReadTestfileInfo($testinfofile);
    $title = $tests{$1}{name};
  }
  if(defined $title)
    {
      FullBox($w);
      TopText($w, $title, CHARSIZE=>1.0, x=>0.5, c=>0.5, y=>-2);
    }
  
}

##########################################

sub SigProp
{
  $PL_nx = 3;
  $PL_ny = 1;
  $PL_ymin = 0.65;
  $PL_ymax = 0.90;
  $PL_deltax = 0.07;
  $cs = 0.5;
  
  $w = NewWindow($psfile);
  
  my $pan = 1;
  
  PlotSigPropCumul($pan++, 3);
  PlotSigPropCumul($pan++, 6);
  PlotSigPropCumul($pan++, 12, 1);
  
  NameTitle($w);  
  $w->close()
}

##########################################

sub Nsig
{
  $w = NewWindow($psfile);
  
  my $pan = 1;
  
  PlotNsig(1);
  $w->close();
}

##########################################

sub PlotNsig
{
  my ($pan) = @_;
  
  my @nsig = ();
  my $maxy = 0;
  my $maxn = 0;
  my @n = (); my @m = (); my @s = ();
  my @med = (); my @med1 = (); my @med2 = ();
  
  local *F;
  open F, $nsigfile or die "Cannot open $nsigfile\n";
  <F>;
  while(my $line = <F>)
  {
    my @d = split /\t/, $line;
    my $n = shift @d;
    my $d = pdl(@d);
    my ($m, $s) = stats($d);
    my ($med, $med1, $med2) = MedianCI($d);
    $nsig[$n] = $d;
    push @n, $n;
    push @m, $m;
    push @s, $s;
    push @med, $med;
    push @med1, $med1;
    push @med2, $med2;
    #print "$med \($med1, $med2)\n";
    my $max = max($d);
    $maxy = $max if $max > $maxy;
    $maxn = $n if $n > $maxn;
  }
    
  $maxy = 1;
  PlotPanelN($w, $pan,
    BOX => [0, $maxn+1, 0, $maxy],
    XLAB => 'Number of replicates',
    YLAB => 'Fraction of significant genes',
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    forceylab => 1,
    YTICK => 0, NYSUB => 0,
    CHARSIZE => $cs
  );
  PlotPanelN($w, $pan,
    BOX => [0, $maxn+1, 0, $maxy],
    XBOX => 'C', YBOX =>'C'
  );
  gridlines($w, 0, $maxn, 5, 0, 1, 0.1);

  my $n = pdl(@n);
  my $med = pdl(@med);
  my $med1 = pdl(@med1);
  my $med2 = pdl(@med2);
  
  #PointPlot($w, $n, $med);
  #errorbar($w, $n, $med, $med-$med1, $med2-$med);

  BoxPlot($w, \@nsig, cloud=>0, solidwhisk=>1, colour=>[215,215,255]);
  #LinePlot($w, $n, pdl(@m), COLOR=>'BLUE', LINEWIDTH=>3);
  #LinePlot($w, $n, pdl(@m)+pdl(@s), COLOR=>'BLUE', LINEWIDTH=>1);
  #LinePlot($w, $n, pdl(@m)-pdl(@s), COLOR=>'BLUE', LINEWIDTH=>1);
}  

##########################################

sub PlotSigProp
{
  my ($pan, $nrep) = @_;
  
  my ($n, $prop) = rcols $sigpropfile, 0, 2;
  my $sel = which($n == $nrep);
  
  my $p = $prop($sel);
  my ($x, $h) = BuildHistogram($p, min=>0, max=>1, renorm=>1);
  
  PlotPanelN($w, $pan,
    BOX => [0, 1, 0, max($h)],
    %PL_empty
  );

  BarPlot($w, $x, $h, colour=>[215,215,255]);
  HistPlot($w, $x, $h);

  PlotPanelN($w, $pan,
    BOX => [0, 1, 0, max($h)],
    XLAB => 'Significance proportion',
    YLAB => 'Frequency',
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    forceylab => 1,
    CHARSIZE => $cs
  );
  $w->close();
  
}

sub PlotSigPropCumul
{
  my ($pan, $nrep, $legend) = @_;
  
  my @fclim = (0, 1, 2);
  my @col = @contrast_colours;
  
  my ($n, $fc, $prop) = rcols $sigpropfile, 0, 2, 3;
  
  my @c = PlotPanelN($w, $pan,
    BOX => [0, 1, 0, 1],
    YLAB => 'DE proportion',
    XLAB => 'Gene fraction',
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    #forceylab => 1,
    CHARSIZE => $cs
  );
  
  for my $i (0 .. @fclim - 1)
  {
    my $fclim = $fclim[$i];
    my $col = $col[$i];
    my $sel = which(($n == $nrep) & (abs($fc) >= $fclim));
    my $p = $prop($sel);
    my $ps = qsort $p;
    my $x = sequence($ps) / (nelem($ps) - 1);
    LinePlot($w, $x, $ps, LINEWIDTH=>3, COLOR=>$col);
  }
  TopText($w, "n = $nrep", y=>1);
  
  my $yl = 0.5;
  LineLegend($w, \@c, [0.6, 0.99, $yl-0.3, $yl], \@col, \@fclim, 'lg#d2#uFC limit', 0.4) if $legend;
}

##########################################

sub PlotFour
{
  my ($pan, $rates) = @_;
  
  my ($nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetTotData($powerfile);
  selrep($nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn);
  my $s = $tp + $tn + $fp + $fn;
  my $maxy = max(append($tp,$tn));
 # $maxy = 100;
  
  my ($tp1, $tn1, $fp1, $fn1, $stp1, $stn1, $sfp1, $sfn1, $fd, $sfd);
  if($rates)
  {
    ($tp1, $stp1) = rerr($tp, $fn, $stp, $sfn);
    ($fp1, $sfp1) = rerr($fp, $tn, $sfp, $stn);
    ($tn1, $stn1) = rerr($tn, $fp, $stn, $sfp);
    ($fn1, $sfn1) = rerr($fn, $tp, $sfn, $stp);
    ($fd, $sfd) = rerr($fp, $tp, $sfp, $stp);
    $maxy = 1;
  }
  else
    {($tp1, $tn1, $fp1, $fn1, $stp1, $stn1, $sfp1, $sfn1) = ($tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn)}
  
  #my $maxy = max(pdl($tp1, $fp1, $tn1, $fn1)->flat()) * 1.05;
  my $maxn=max($nrep);
  my @c = PlotPanelN($w, $pan,
    BOX => [0, $maxn+1, 0, $maxy],
    %PL_empty
  );
  
  gridlines($w, 0, $maxn+1, 5, 0, $maxy, 500);
  
  my @col = (
    [0, 0, 0],
    [255, 0, 0],
    [0, 127, 0],
    [191, 0, 191]
  );
  
  shadeerr($w, $nrep, $tp1, $stp1, $col[0]);
  line($w, $nrep, $tp1, $col[0]);
  shadeerr($w, $nrep, $fp1, $sfp1, $col[1]);
  line($w, $nrep, $fp1, $col[1]);
  shadeerr($w, $nrep, $tn1, $stn1, $col[2]);
  line($w, $nrep, $tn1, $col[2]);
  shadeerr($w, $nrep, $fn1, $sfn1, $col[3]);
  line($w, $nrep, $fn1, $col[3]);
  
  LineLegend($w, \@c, [0.8, 0.95, 0.6, 0.8], \@col, ['TP', 'FP', 'TN', 'FN'], undef, $cs2);
  
  PlotPanelN($w, $pan,
    BOX => [0, $maxn+1, 0, $maxy],
    XLAB => 'Number of replicates',
    YLAB => '',
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    forceylab => 1,
    YTICK => 0, NYSUB => 0,
    CHARSIZE => $cs
  );
  PlotPanelN($w, $pan,
    BOX => [0, $maxn+1, 0, $maxy],
    XBOX => 'C', YBOX =>'C'
  );
  
  
  #linewitherr($w, $nrep, $tp1, $stp1, 'BLACK');
  #linewitherr($w, $nrep, $fp1, $sfp1, 'RED');
  #linewitherr($w, $nrep, $tn1, $stn1, 'BLUE');
  #linewitherr($w, $nrep, $fn1, $sfn1, 'GOLD2');
  #linewitherr($w, $nrep, $fd, $sfd, 'GREEN') if $rates;
  
  my $tcs = 0.4;
  my $rt = ($rates) ? 'R' : '';
  my ($xx, $dx) = (0.1, 0.15);
  #TopText($w, "TP$rt", y=>1, x=>$xx, COLOR=>$col[0], CHARSIZE=>$tcs);
  #TopText($w, "FP$rt", y=>1, x=>$xx+$dx, COLOR=>$col[1], CHARSIZE=>$tcs);
  #TopText($w, "TN$rt", y=>1, x=>$xx+2*$dx, COLOR=>$col[2], CHARSIZE=>$tcs);
  #TopText($w, "FN$rt", y=>1, x=>$xx+3*$dx, COLOR=>$col[3], CHARSIZE=>$tcs);
}
  
##########################################

sub PlotFC
{ 
  my ($pan, $what) = @_;
  
  my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetFCData($fcfile1);
  selrep($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn);
  
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

  PlotPanelN($w, $pan,
    BOX => [-4, 3, 0, 1],
    XLAB => 'log2 FC',
    YLAB => $lab,
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    forceylab => 1,
    YTICK => 0, NYSUB => 0,
    CHARSIZE => $cs
  );

  
  my $i = 0;
  my $x = -3.8;
  for my $n (@nrep)
  {
    my @fc = (); my @R = (); my @sR = ();
    for my $f (@F)
    {
      my $sel = which(($fclim == $f) & ($nrep == $n));
      next if nelem($sel) == 0;
      #print "S = $sel\n";
      my ($rf, $srf) = ($r($sel), $sr($sel));
      push @fc, $f;
      push @R, at($rf, 0);
      push @sR, at($srf, 0);
    }
    my $col = ($i >= @col) ? 'GREY' : $col[$i];
    #linewitherr($w, pdl(@fc), pdl(@R), pdl(@sR), $col);
    LinePlot($w, pdl(@fc), pdl(@R), COLOR=>$col, LINEWIDTH=>3);
    #wcols pdl(@fc), pdl(@R), pdl(@sR);
    my $y = 0.88 - 0.06 * $i;
    LinePlot($w, pdl($x,$x+0.2), pdl($y,$y), COLOR=>$col, LINEWIDTH=>3);
    $w->text($n, COLOR=>$col, CHARSIZE=>$cs2, TEXTPOSITION=>[$x+0.3,$y,0,0,0]);
    
    $i++
  }
  $w->text('N rep', COLOR=>'BLACK', CHARSIZE=>$cs2, TEXTPOSITION=>[$x,0.96,0,0,0]);
}


##########################################

sub PlotNrepRateFC
{ 
  my ($pan, $what, %opt) = @_;
  
  my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetFCData($fcfile2);
  selrep($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn);
  my $maxn = max($nrep);
  my $maxy = 1;

  my ($lab, $r, $sr);
  if($what eq 'fnr')
  {
    $lab = 'FNR';
    ($r, $sr) = rerr($fn, $tp, $sfn, $stp)
  }
  elsif($what eq 'fpr')
  {
    $lab = 'False positive rate';
    ($r, $sr) = rerr($fp, $tn, $sfp, $stn);
    $maxy = 0.2;
  }
  elsif($what eq 'tpr')
  {
    $lab = 'True positive rate';
    ($r, $sr) = rerr($tp, $fn, $stp, $sfn)
  }
  #elsif($what eq 'sigfrac')
  #{
   # $lab = "Significant fraction";
    #($r, $sr) = ($tp / $nsig, $stp / $nsig)
  #}
  elsif($what eq 'fdr')
  {
    $lab = 'FDR';
    ($r, $sr) = rerr($fp, $tp, $sfp, $stp);
    $maxy = 0.2;
  }
  $lab = "True (false) positive rate" if $opt{withfalse};
  

#wcols $nrep, $fclim, $fn, $tp, $r, $sr;die;$nsig + $nnonsig

  my $miny = (defined $opt{miny}) ? $opt{miny} : 0;
  $maxy = $opt{maxy} if defined $opt{maxy};
  $maxn = $opt{maxn} if defined $opt{maxn};

  my @c = PlotPanelN($w, $pan,
    BOX => [0, $maxn+1, $miny, $maxy], %PL_empty
  );
  
  my $gy = ($maxy < 0.4) ? 0.05 : 0.1;
  gridlines($w, 0, $maxn, 5, 0, 1, $gy);

  
  my @x = ();
  my @y = ();
  my @sy = ();
  
  my $i = 0;
  my $x = $maxn+2;
  for my $f (@fclim)
  {
    my @n = (); my @R = (); my @sR = ();
    for my $n (2 .. $maxn)
    {
      my $sel = which(($fclim == $f) & ($nrep == $n));
      next if nelem($sel) == 0;
      #print "S = $sel\n";
      my ($rf, $srf) = ($r($sel), $sr($sel));
      push @n, $n;
      push @R, sum($rf);
      push @sR, sqrt(sum($srf**2));
    }
    my $col = ($i >= @col) ? 'GREY' : $col[$i];
    #linewitherr($w, pdl(@n), pdl(@R), pdl(@sR), $col);
    push @x, pdl(@n);
    push @y, pdl(@R);
    push @sy, pdl(@sR);
    #LinePlot($w, pdl(@n), pdl(@R), COLOR=>$col, LINEWIDTH=>3);
    #wcols pdl(@n), pdl(@R), pdl(@sR);die;
    #my $y = 0.88 - 0.06 * $i;
    #LinePlot($w, pdl($x,$x+0.9), pdl($y,$y), COLOR=>$col, LINEWIDTH=>3);
    #$w->text($f, COLOR=>$col, CHARSIZE=>$cs2, TEXTPOSITION=>[$x+1.5,$y,0,0,0]);
    #$nsig + $nnonsig
    $i++
  }
  #$w->text('lg#d2#uFC', COLOR=>'BLACK', CHARSIZE=>$cs2, TEXTPOSITION=>[$x,0.96,0,0,0]);
  
  shadeerr($w, $x[$_], $y[$_], $sy[$_], $col[$_]) for (0 .. @x - 1);
  line($w, $x[$_], $y[$_], $col[$_]) for (0 .. @x - 1);
  
  if($opt{withfalse})
  {
    my ($nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetTotData($powerfile);
    selrep($nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn);
    my $maxn = max($nrep);

    my ($fpr, $sfpr) = rerr($fp, $tn, $sfp, $stn);
    my $fcol = $col[0];
  
    shadeerr($w, $nrep, $fpr, $sfpr, $fcol);
    line($w, $nrep, $fpr, $fcol, 4);
  }
    
  my $yl = (defined $opt{ylegend}) ? $opt{ylegend} : 0.5;
  LineLegend($w, \@c, [0.7, 0.92, $yl-0.3, $yl], \@col, \@fclim, '|lg#d2#uFC|>', $cs2);
  
  PlotPanelN($w, $pan,
    BOX => [0, $maxn+1, $miny, $maxy],
    XLAB => 'Number of replicates',
    YLAB => $lab,
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    forceylab => 1,
    YTICK => 0, NYSUB => 0,
    CHARSIZE => $cs,
  );
  PlotPanelN($w, $pan,
    BOX => [0, $maxn+1, $miny, $maxy],
    XBOX => 'C', YBOX =>'C'
  );
  
  
}

##########################################

sub PlotNrepRate2
{ 
  my ($pan, %opt) = @_;
  
  my ($nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetTotData($powerfile);
  selrep($nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn);
  my $maxn = max($nrep);

  my ($fpr, $sfpr) = rerr($fp, $tn, $sfp, $stn);
  my ($tpr, $stpr) = rerr($tp, $fn, $stp, $sfn);

  my ($miny, $maxy) = (0, 1);

  my @c = PlotPanelN($w, $pan,
    BOX => [0, $maxn, $miny, $maxy],
    XLAB => 'Number of replicates',
    YLAB => 'True/false positive rate',
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    forceylab => 1,
    YTICK => 0, NYSUB => 0,
    CHARSIZE => $cs,
  );
  
  gridlines($w, 0, $maxn, 5, 0, 1, 0.1);
  
  my ($tcol, $fcol) = ([0,0,255], [255,0,0]);
  
  shadeerr($w, $nrep, $tpr, $stpr, $tcol);
  line($w, $nrep, $tpr, $tcol);
  shadeerr($w, $nrep, $fpr, $sfpr, $fcol);
  line($w, $nrep, $fpr, $fcol);
  
  LineLegend($w, \@c, [0.7, 0.9, 0.3, 0.5], [$tcol,$fcol], ['TPR','FPR'], undef, 0.4);
}

##########################################

sub PlotNrepFrac
{ 
  my ($pan, $what) = @_;
  
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

#wcols $nrep, $fclim, $fn, $tp, $r, $sr;die;

  PlotPanelN($w, $pan,
    BOX => [0, $maxn+7, 0, 1],
    XLAB => 'Number of replicates',
    YLAB => $lab,
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    forceylab => 1,
    YTICK => 0, NYSUB => 0,
    CHARSIZE => $cs
  );

  
  my $i = 0;
  my $x = $maxn+2;
  for my $f (@fclim)
  {
    my @n = (); my @R = (); my @sR = ();
    for my $n (2 .. $maxn)
    {
      my $sel = which(($fclim == $f) & ($nrep == $n));
      next if nelem($sel) == 0;
      #print "S = $sel\n";
      my ($rf, $srf) = ($r($sel), $sr($sel));
      push @n, $n;
      push @R, sum($rf);
      push @sR, sqrt(sum($srf**2));
    }
    my $col = ($i >= @col) ? 'GREY' : $col[$i];
    linewitherr($w, pdl(@n), pdl(@R), pdl(@sR), $col);
    #LinePlot($w, pdl(@n), pdl(@R), COLOR=>$col, LINEWIDTH=>3);
    #wcols pdl(@n), pdl(@R), pdl(@sR);die;
    my $y = 0.88 - 0.06 * $i;
    LinePlot($w, pdl($x,$x+0.9), pdl($y,$y), COLOR=>$col, LINEWIDTH=>3);
    $w->text($f, COLOR=>$col, CHARSIZE=>$cs2, TEXTPOSITION=>[$x+1.5,$y,0,0,0]);
    
    $i++
  }
  $w->text('lg#d2#uFC', COLOR=>'BLACK', CHARSIZE=>$cs2, TEXTPOSITION=>[$x,0.96,0,0,0]);
}

##########################################

sub PlotFCRateNrep
{ 
  my ($pan, $what, %opt) = @_;
  
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
    $lab = 'False positive rate';
    ($r, $sr) = rerr($fp, $tn, $sfp, $stn)
  }
  elsif($what eq 'tpr')
  {
    $lab = 'True positive rate';
    ($r, $sr) = rerr($tp, $fn, $stp, $sfn)
  }
  #elsif($what eq 'sigfrac')
  #{
  #  $lab = "Significant fraction";
   # ($r, $sr) = ($tp / $nsig, $stp / $nsig)
 # }
  $lab = "True (false) positive rate" if $opt{withfalse};
  

  #wcols $nrep, $fclim, $fn, $tp, $r, $sr;die;

  my @c = PlotPanelN($w, $pan,
    BOX => [0, $maxfc, 0, 1], %PL_empty
  );
  gridlines($w, 0, 3, 0.5, 0, 1, 0.1);

  my @x = ();
  my @y = ();
  my @sy = ();

  
  my $i = 0;
  my $x = 2.5;
  for my $n (@nrep)
  {
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
    my $col = ($i >= @col) ? 'GREY' : $col[$i];
    #linewitherr($w, pdl(@f), pdl(@R), pdl(@sR), $col);
    push @x, pdl(@f);
    push @y, pdl(@R);
    push @sy, pdl(@sR);
    
    #LinePlot($w, pdl(@n), pdl(@R), COLOR=>$col, LINEWIDTH=>3);
    #wcols pdl(@n), pdl(@R), pdl(@sR);die;
    #my $y = 0.5 - 0.06 * $i;
    #LinePlot($w, pdl($x,$x+0.2), pdl($y,$y), COLOR=>$col, LINEWIDTH=>3);
    #$w->text($n, COLOR=>$col, CHARSIZE=>$cs2, TEXTPOSITION=>[$x+.3,$y,0,0,0]);
    
    $i++
  }
  #$w->text('N#drep', COLOR=>'BLACK', CHARSIZE=>$cs2, TEXTPOSITION=>[$x,0.6,0,0,0]);
  
  shadeerr($w, $x[$_], $y[$_], $sy[$_], $col[$_]) for (0 .. @x - 1);
  line($w, $x[$_], $y[$_], $col[$_]) for (0 .. @x - 1);
  
  if($opt{withfalse})
  {
    my ($r, $sr) = rerr($fp, $tn, $sfp, $stn);
    my @f = (); my @R = (); my @sR = ();
    for my $f (@F)
    {
      my $sel = which(($fclim == $f) & ($nrep == 3));
    #print "FP=", $fp($sel), "\nTN=", $tn($sel), "\n\n";
      next if nelem($sel) == 0;
      my ($rf, $srf) = ($r($sel), $sr($sel));
      my $sumrf = sum($rf);
      unless($sumrf =~ /nan/)
      {
        push @f, $f;
        push @R, sum($rf);
        push @sR, sqrt(sum($srf**2));
        #print $f, " ", $rf, "\n" unless sum($rf) =~ /nan/;;
      }
    }
    my $fcol = $col[0];
    shadeerr($w, pdl(@f), pdl(@R), pdl(@sR), $fcol);
    line($w, pdl(@f), pdl(@R), $fcol, 4);
  }
  
  LineLegend($w, \@c, [0.7, 0.9, 0.15, 0.5], \@col, \@nrep, 'N#drep#u', $cs2);
  
  PlotPanelN($w, $pan,
    BOX => [0, $maxfc, 0, 1],
    XLAB => 'log#d2#uFC',
    YLAB => $lab,
    xbox => 'NI', ybox => 'NI',
    forcexlab => 1,
    forceylab => 1,
    YTICK => 0, NYSUB => 0,
    CHARSIZE => $cs
  );
  PlotPanelN($w, $pan,
    BOX => [0, $maxfc, 0, 1],
    XBOX => 'C', YBOX =>'C'
  );
  
  
}

#####################################

sub PlotROC
{
  my ($pan) = @_;
  
  my ($nrep, $fclim, $fpr, $tpr, $sfpr, $stpr) = rcols $rocfile, {LINES=>'1:-1'};
  my $sel = which(($nrep == $nroc) & ($fclim == $froc));
  $fpr = $fpr($sel);
  $tpr = $tpr($sel);
  $sfpr = $sfpr($sel);
  $stpr = $stpr($sel);
  
  PlotPanelN($w, $pan,
    BOX => [0, 1, 0, 1],
    XLAB => 'FPR',
    YLAB => 'TPR',
    forcexlab=>1, forceylab=>1,
    XBOX=>'BNSTI', YBOX=>'BNSTI',
    CHARSIZE=>$cs
  );
  
  LinePlot($w, pdl(0,1), pdl(0,1), COLOR=>'RED');
  LinePlot($w, $fpr, $tpr, COLOR=>'GREY');
  PointPlot($w, $fpr, $tpr, XERRORBAR=>2*$sfpr, YERRORBAR=>2*$stpr, COLOR=>'GREY', SYMBOLSIZE=>0.001, MINTICKSIZE=>0.001);
  PointPlot($w, $fpr, $tpr);
  TopText($w, "|lg#d2#uFC| > $froc, n#drep#u = $nroc", CHARSIZE=>0.4, y=>1);
  
}
  
sub line
{
  my ($w, $x, $r, $col, $sty) = @_;
  
  $sty ||= 1;
  LinePlot($w, $x, $r, LINEWIDTH=>3, COLOR=>$col, LINESTYLE=>$sty);
}


sub linewitherr
{
  my ($w, $x, $r, $sr, $col) = @_;
  LinePlot($w, $x, $r, LINEWIDTH=>3, COLOR=>$col);
  LinePlot($w, $x, $r+$sr, LINEWIDTH=>1, COLOR=>$col);
  LinePlot($w, $x, $r-$sr, LINEWIDTH=>1, COLOR=>$col);
}

sub shadeerr
{
  my ($w, $x, $r, $sr, $col) = @_;
  
  my $px = append($x, $x(-1:0));
  my $py = append($r - $sr, $r(-1:0) + $sr(-1:0));
  
  my $shade;
  my $fade = 0.2;
  $shade->[$_] = 255 - (255 - $col->[$_]) * $fade for (0 .. 2);
  
  plscmap1(pdl($shade->[0],0), pdl($shade->[1],0), pdl($shade->[2],0));
  plcol1(0);
  plfill($px, $py);

 # LinePlot($w, $x, $r, LINEWIDTH=>3, COLOR=>$col);
}

sub gridlines
{
  my ($w, $x1, $x2, $xs, $y1, $y2, $ys) = @_;
  
  my $col = 'GREY';
  for(my $x = $x1; $x <= $x2; $x += $xs)
    {LinePlot($w, pdl($x,$x), pdl(-1e4,1e4), COLOR=>$col)}
  for(my $y = $y1; $y <= $y2; $y += $ys)
    {LinePlot($w, pdl(-1e4,1e4), pdl($y,$y), COLOR=>$col)}
}


sub LineLegend
{
  my ($w, $wcoord, $lcoord, $cols, $labs, $lab, $cs) = @_;
  $cs ||= 0.5;
  
  my ($x1, $x2, $y1, $y2) = @$wcoord;
  my ($fx1, $fx2, $fy1, $fy2) = @$lcoord;
  
  my $wx1 = $x1 + $fx1 * ($x2 - $x1);
  my $wx2 = $x1 + $fx2 * ($x2 - $x1);
  my $wy1 = $y1 + $fy1 * ($y2 - $y1);
  my $wy2 = $y1 + $fy2 * ($y2 - $y1);
  
  $w->xyplot(pdl(-1),pdl(-1),
    VIEWPORT   => [$wx1, $wx2, $wy1, $wy2],
    BOX => [0, 1, 0, 1],
    %PL_empty
  );
  
  my $px = pdl(0,1,1,0,0);
  my $py = pdl(0,0,1,1,0);
  plscmap1(pdl(255,0), pdl(255,0), pdl(255,0));
  plcol1(0);
  plfill($px, $py);
  LinePlot($w, $px, $py, COLOR=>'GREY');
  
  my $n = @$labs + 1;
  my $step = (defined $lab) ? 1 / ($n + 0.5) : 1 / ($n - 1);
  my $left = 0.1;
  my $div1 = 0.5;
  my $div2 = 0.6; 
  
  $w->text($lab, COLOR=>'BLACK', CHARSIZE=>$cs, TEXTPOSITION=>[$left,1 - $step/2,0,0,0]) if defined $lab;
  
  my $di = (defined $lab) ? 1.7 : 0.5;
  for my $i (0 .. $n - 1)
  {
    my $l = $labs->[$i];
    my $c = $cols->[$i];
    $c ||= 'GREY';
    
    my $y = 1 - ($i + $di) * $step;
    LinePlot($w, pdl($left, $div1), pdl($y, $y), COLOR=>$c, LINEWIDTH=>3);
    $w->text($l, COLOR=>$c, CHARSIZE=>$cs, TEXTPOSITION=>[$div2,$y,0,0,0]);
  }
  
}


sub selrep
{
  my $nrep = $_[0];
  my $sel = which($nrep >= $minnrep);
  for my $i (0 .. @_ - 1)
    {$_[$i] = $_[$i]->($sel)}
}

=head1 SYNOPSIS

  plot_powerstats_db.pl - test=edger
    
=head1 OPTIONS

=over 4

=item B<-graph>=[over|sim|cn|ns|sp]

Which graph to plot. Only 'over' and 'sp' are working.

  over - overview plot, used in the paper (default)
  sp - significance propotion for three replicate numbers


=item B<-testinfofile>=I<pathname>

A file with DE tools metadata. See C<plot_fcp.pl> for details. Default name is C<de_tests.txt>

=item B<-powerdir>=I<string>

Directory where the results from C<make_powerstats_db.pl> are stored. The default value is C<powerstats_db_ref>.


=item B<-test>=I<string>

Name of the test/tool to plot results for. As defined in the metadata file.


=item B<-minnerp>=I<integer>

Minumum number of replicates to include in the plots. Default is 2.


=item B<-maxfc>=I<number>

Maximum value of log2-fold-change axis in the plots. The default value is 2.  

=item B<-cs>=I<number>

Larger character size in plots. The default value is 0.45. 

=item B<-cs2>=I<number>

Smaller character size in plots. The default value is 0.3.

=item B<-psfile>=I<pathname>

Name of the output postscript file.


=back
