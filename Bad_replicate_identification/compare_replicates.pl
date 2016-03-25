#!/sw/bin/perl

=head1 NAME

B<compare_replicates.pl>

=head1 DESCRIPTION
  
Show various comparisons between replicates in a given condition.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use GRNASeq;

use Stats;
use PLgraphs;
use Tools;
use Heatmap;

use Switch 'Perl6';

ReadDefs();
Stamp();

$| = 1;

my $psfile;
my $cond = 'WT';
my ($rep, $rep1, $rep2);
my $graph = 'all';
my $norm = 'deseq';
my $type = 'raw';
my $log;
my $exclude;
my $outliers;
my $genlist;
my $nonzero;
my ($help, $man);
GetOptions(
  'psfile=s' => \$psfile,
  #'rep=i' => \$rep,
  'rep1=s' => \$rep1,
  'rep2=s' => \$rep2,
  'cond=s' => \$cond,
  'graph=s' => \$graph,
  'norm=s' => \$norm,
  #'exclude=s' => \$exclude,
  #'outliers=i' => \$outliers,
  #'genlist=s' => \$genlist,
  #'nonzero=i' => \$nonzero,
  #'type=s' => \$type,
  #log => \$log,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;


my %reps = ( 
  WT => ['4,27,7,4,10,10,10,25,21', '48,31,33,41,43,22,28,34,34'],
  Snf2 => ['5,5,16,15,14,21,3,6,6', '43,26,48,27,25,35,13,42,10']
);

my $cs = 0.5;

my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel, $f);
my %g2i;
my $N;

unless($graph eq 'all' || $graph eq 'mc')
{
  ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel, $f) = GetNormalizedData(
    condition=>$cond,
    replicate=>$rep,
    exclude=>$exclude,
    nonzero=>$nonzero,
    norm=>$norm,
    type=>$type
  );
  %g2i = GenHash($genes, $sel);
  $N = sumover($d>0);
}

if(defined $outliers)
{
  my $rr = sprintf "%02d", $outliers;
  my $limit = 5;
  $genlist = "$outlierdir/${cond}_${norm}_rep${rr}_s${limit}.dat";
  $name .= " with outliers rep$rr";
}

my $glist = GeneIdList($genlist, \%g2i) if defined $genlist;

my @cl = split(/:/, $CleanExclude);
my %clean = ();
for my $i (0 .. 1)
{
  my @c = split(/,/, $cl[$i]);
  $clean{$conds[$i]}{$_} = 1 for (@c);
}


given($graph)
{
  when 'all'   {AllReplicates()}
  when '2'   {TwoReplicates()}
  #when 'cd' {CorrDist()}
  when 'rm' {RepMeanMany()}
  when 'mc' {MeanCorr()}
}

##############################################



sub AllReplicates
{
  my $w = NewWindow($psfile);
    
  my $cs = 0.7;
    
  $PL_ymin = 0.50;
  $PL_ymax = 0.95;
  $PL_xmax = 0.85;
  $PL_deltax = 0.02;
  $PL_deltay = 0;
  $PL_nx = 2;
  my ($min, $max) = (0.8, 1);  

  my $C;
  
  my $pan = 1;
  for my $cond (@conds)
  {
    my ($d, $m, $s, $nrep, $ngen) = GetNormalizedData(
      condition=>$cond,
      nonzero=>1,
      norm=>'none',
      type=>$type
    );
    $C = CorrelationMatrix($d);
    
    PlotPanelN($w, $pan,
      BOX=>[0, $nrep-1, 0, $nrep-1],
      %PL_empty   
    );
    
    SCMap('heat');
    HeatMapRows($C, {min=>$min, max=>$max});
    PlotPanelN($w, $pan,
      BOX=>[0.5, $nrep+0.5, 0.5, $nrep+0.5],
      xbox=>'CI', ybox=>'CI',
      XLAB => 'Replicate', YLAB => 'Replicate',
      CHARSIZE => $cs
    );
    TopText($w, $cond, x=>0, y=>1);
    
    $pan++;
  }
  
  $w->colorkey($C, 'v', 
    VIEWPORT => [0.87, 0.90, $PL_ymin, $PL_ymax],
    YBOX => 'CMSTI',
    XLAB => '', YLAB => '',
    #YLAB => 'Correlation',
    CHARSIZE => $cs,
    ZRANGE => [$min, $max]
  );
  
  #FullBox($w);
  #TopText($w, $name, x=>0.07, y=>-2, CHARSIZE => 0.7);
  #TopText($w, "Normalization: $norm", x=>0.07, y=>-4, CHARSIZE => 0.7);
  
  $w->close();
}

sub MeanCorr
{
  my $w = NewWindow($psfile);
    
  my $cs = 0.7;
    
  $PL_ymin = 0.60;
  $PL_ymax = 0.95;
  $PL_xmax = 0.85;
  $PL_deltax = 0.02;
  $PL_deltay = 0;
  $PL_nx = 2;
  my ($min, $max) = (0.8, 1);  
  
  my $pan = 1;
  for my $cond (@conds)
  {
    my ($d, $m, $s, $nrep, $ngen) = GetNormalizedData(
      condition=>$cond,
      nonzero=>1,
      norm=>'none',
      type=>$type
    );
    my $C = CorrelationMatrix($d);
    my ($MC, $MS, $MMED) = statsover($C); 
    
    PlotPanelN($w, $pan,
      BOX=>[0, $nrep+1, 0.65, 1],
      xbox=>'I', ybox=>'I',
      XLAB => 'Replicate', YLAB => 'Median correlation',
      CHARSIZE => $cs
    );
    
    PointPlot($w, sequence($MMED)+1, $MMED, SYMBOLSIZE=>1);
    LinePlot($w, pdl($_,$_), pdl(1,at($MMED, $_-1)), COLOR=>'GREY') for (1..$nrep);
    
    for my $i (0 .. nelem($MC)-1)
    {
      my $mc = at($MMED, $i);
      if($clean{$cond}{$i+1})
        {$w->text($i+1, CHARSIZE=>0.4, COLOR=>'RED', TEXTPOSITION=>[$i+1, $mc-0.01, 0, 0, 0.5])}
    }
    
    TopText($w, $cond, x=>0, y=>1);
    
    $pan++;
  }
  
  #FullBox($w);
  #TopText($w, $name, x=>0.07, y=>-2, CHARSIZE => 0.7);
  #TopText($w, "Normalization: $norm", x=>0.07, y=>-4, CHARSIZE => 0.7);
  
  $w->close();
  
}

sub CorrDist
{
  my $R = CorrelationMatrix();;
  #$R = $R->where($R>0.97);
  my ($x, $h) = BuildHistogram($R, renorm=>1, extend=>1);
  
  my $w = NewWindow($psfile);
  PlotPanelN($w, 1,
    BOX => [0.5, 1.0, 0, max($h)*1.02],
    XLAB => 'Correlation coefficient',
    YLAB => 'Normalized frequency'
  );
  
  BarPlot($w, $x, $h, colour=>[215,215,255]);
  HistPlot($w, $x, $h);
  
  $w->close();  
  
}


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
      my $r = ($log) ? corr(log10($r1($sel)), log($r2($sel))) 
                     : corr($r1($sel), $r2($sel));
                     
    #  $r = logrms($r1, $r2);
                     
      #printf "%02d %02d %8.6f\n", $i+1, $j+1, $r;
      $C($i, $j) .= $r;
      $C($j, $i) .= $r;
    }
  }
  return $C;
}

sub TwoReplicates
{
  
  ($rep1, $rep2) = @{$reps{$cond}} unless(defined $rep1 && defined $rep2);
  
  my @r1 = split(/,/, $rep1);
  my @r2 = split(/,/, $rep2);
  die "Need equal number of replicates\n" unless scalar @r1 == scalar @r2;
  my $n = @r1;
  
  my $nx = ceil(sqrt($n));
  $PL_nx = $nx;
  $PL_ny = $nx;
  $PL_deltax = 0.05;
  $PL_deltay = 0.05;
  $PL_xmin = 0.1;
  $PL_ymin = 0.1;
  $cs = 1.1;
  
  my $w = NewWindow($psfile);
  
  for my $i (0 .. @r1 - 1)
    {TwoRepPanel($w, $i+1, $r1[$i], $r2[$i])}  
  
  FullBox($w);
  TopText($w, "Replicate comparison for $cond ($norm)", x=>0.5, c=>0.5);
  
  $w->close();
}


sub TwoRepPanel
{
  my ($w, $pan, $rr1, $rr2) = @_;
  
  my $r1 = $d($rr1-1,;-);
  my $r2 = $d($rr2-1,;-);
  
  my $sel = which(($r1 > 0) & ($r2 > 0));
  #my $r = corr($r1($sel), $r2($sel));
  
  $r1 = $r1($sel);
  $r2 = $r2($sel);
  
  
  $r1->where($r1==0) .= 0.1;
  $r2->where($r2==0) .= 0.1;
  
  my $l1 = log10($r1);
  my $l2 = log10($r2);
  
  my $max = max($l1->append($l2));
  
#  PlotPanelN($w, $pan,
#    BOX => [-1.2, $max, -1.2, $max],
#    XLAB => "log Rep1",
#    YLAB => "log Rep2",
#    CHARSIZE => $cs,
#  );
#  
#  my $p = which(($r1 > 0.1) & ($r2 > 0.1));
#  my $n = which(($r1 == 0.1) | ($r2 == 0.1));
#    
#  my $col = ($genlist) ? 'GREY' : 'BLACK';
#  PointPlot($w, $l1($p), $l2($p), SYMBOLSIZE=>0.5, COLOR=>$col);
#  PointPlot($w, $l1($glist), $l2($glist), SYMBOLSIZE=>0.5, COLOR=>'BLACK') if defined $glist;
#  PointPlot($w, $l1($n), $l2($n), SYMBOLSIZE=>0.5, COLOR=>'GOLD2');
#  LinePlot($w, pdl(-3,$max), pdl(-3,$max), COLOR=>'RED');

  my ($F, $r) = MAPlot($w, $pan, $r1, $r2, ymin=>-6, ymax=>6, xmin=>-1, xmax=>6, nbin=>300, plottype=>'shade');
  TopText($w, sprintf("#gF=%.2g", $F), x=>1, c=>1);
  TopText($w, sprintf("r=%.2g", $r), x=>1, c=>1, y=>-3);
  
  TopText($w, "$rr1:$rr2", CHARSIZE=>0.8, y=>-1);
  #TopText($w, sprintf("r = %6.4f", $r), CHARSIZE=>0.5, y=>-3.5);
  
}



sub RepMeanMany
{
  $PL_nx = 4;
  $PL_ny = 4;
  
  for my $page (1 .. 3)
  {
    my $w = NewWindow("${cond}_${norm}_meanrep_$page.ps");  
    my $pan = 1;
    my $r1 = ($page - 1) * $PL_nx * $PL_ny + 1;
    my $r2 = $r1 + $PL_nx * $PL_ny - 1;
    for my $r ($r1 .. $r2)
    {
      #print "$r ";
      RepMean($w, $pan++, $r); 
    }  
    
    FullBox($w);
    $name .= " norm: $norm";
    TopText($w, $name, x=>0.07, y=>-2, CHARSIZE => 0.7);
    $w->close();
  }
  print "\n";
}

sub RepMean
{
  my ($w, $pan, $rep) = @_;
  
  my $r = $d($rep-1,;-);
  my ($m, $s) = statsover($d);
  
  my $sel = which(($r > 0) & ($m > 0));
  $r = $r($sel);
  $m = $m($sel);
  my $ratio = $r / $m;
  my ($mr, $sr) = stats($ratio);
  my $vr = $sr**2;
  print "$rep: $vr\n";
  #print $r->where($r<10), "\n";
  
  my $lm = log10($m);
  my $lr = log($ratio) / log(2);
  

  PlotPanelN($w, $pan,
    BOX => [-1.5, 6, -3.999, 4],
    XLAB => "log M",
    YLAB => "log#d2#u (R/M)",
  );
  
  PointPlot($w, $lm, $lr, SYMBOLSIZE=>0.1);
  LinePlot($w, pdl(-2,6), pdl(0,0), COLOR=>'RED');
  TopText($w, $rep, c=>1, x=>1, CHARSIZE=>0.8);
}
 
 
sub logrms
{
  my ($x, $y) = @_;
  
  my $sel = which(($x > 0) & ($y > 0));
  my $rat = log($y($sel) / $x($sel)) / log(2);
  
  return sqrt(avg($rat**2));
}



=head1 SYNOPSIS

  compare_replicates.pl -graph=all
  
=head1 OPTIONS

There are several other options in this script, used for exploratory analysis, not actually used in the paper.

=over 4

=item B<-graph>=I<[all|2]>
 
Type of plot to create.

  * all - heat plots with all replicates
  * 2 - two replicates (or selection of pairs)

=item B<-rep1>=I<number>

First replicate to compare

=item B<-rep2>=I<number>

Second replicate to compare

=item B<-cond>=I<string>

Condition name.

=item B<-cond>=I<string>

Condition name (i.e. WT or Snf2).

=item B<-psfile>=I<pathname>

Output postscript file.


=back

