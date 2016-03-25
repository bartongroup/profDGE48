#!/sw/bin/perl

=head1 NAME

B<poisson_reps.pl>

=head1 DESCRIPTION

Performs goodness-of-fit test for Poisson distribution across lanes. Requires gene expression files in a default directory (as defined in C<defs.dat>) and a total read count file created with C<count_fastq_reads_lane.pl>.

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
use Distribution;
use PLgraphs;
use Tools;

$| = 1;

Stamp();
ReadDefs();


my %distribution = (
  norm => 'lognormal',
  pois => 'Poisson'
);

my %testtype = (
  dispersion => 'dispersion',
  logratio => 'log ratio'
);

my $readcountfile = "all_readcount_lane.dat";
my $outfile = 'all_poiss.dat';
my $psfile;
my $type = 'lev';
my $ptype = 'dispersion';
my $alpha = 0.05;
my $myminp;
my ($help, $man);
GetOptions(
  'readcountfile=s' => \$readcountfile,
  'outfile=s' => \$outfile,
  #'type=s' => \$type,
  'ptype=s' => \$ptype,
  #'minp=f' => \$myminp,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;


my ($minp, $maxp) = (-17.1, 0);
$minp = $myminp if defined $myminp; 

my %lanetot = GetLaneTot() unless $type eq 'lev';

my $tot = 0;
my $sig = 0;

my ($P, $M, $S, $T);
my $first = 1;

for my $cond (@conds)
{
  for my $rep (1 .. 48)  
  {
    my $d = GetRep($cond, $rep);
    my ($nlan, $ngen) = dims $d;
    my $p = PoissonFromData($d, $ptype, 1);
    my ($m, $s) = statsover($d);
    my $t = ($nlan - 1) * $s**2 / $m;        # T statistic, chi-square dist
    
    if($first)
    {
      $P = $p;
      $M = $m;
      $S = $s;
      $T = $t;
      $first = 0;
    }
    else
    {
      $P = append($P, $p);
      $M = append($M, $m);
      $S = append($S, $s);
      $T = append($T, $t);
    }
    
    my $n = nelem($p);
    #my $limit = HolmBonferroniLimit($P, $alpha);
    my $limit = $alpha / $n; 
    
    my $ns = nelem($p->where($p < $limit));
    my $f = $ns / $n;
    printf "rep = %d  limit = %.3g  Ns = %d  f = %.3g\n", $rep, $limit, $ns, $f;
    
    $tot += $n;
    $sig += $ns;
    
    #for my $i (0 .. $n-1)
    #  {print $d(,$i;-), "   ", at($P, $i), "\n"}
    #die;
    
    #PlotPDist($w, $rep, $P, $ns, $limit);
  }
}

wcols $M, $S, $P, $T, $outfile;
print "Results written in $outfile\n";

my $ntot = nelem($P);
my $limit = $alpha / $ntot; 
my $ns = nelem($P->where($P < $limit));

my ($hblimit, $nshb) = HolmBonferroniLimit($P, $alpha);

print "limit = $limit\n";
print "Total = $ntot, significant = $ns  hb = $nshb\n";
print "Total ind = $tot, significant ind = $sig\n";



#################################


sub GetLaneTot
{
  my %h = ();
  open F, $readcountfile or die "Cannot open file $readcountfile\n";
  while(my $line = <F>)
  {
    chomp $line;
    my ($cond, $rep, $lane, $f, $count) = split " ", $line;
    $h{$cond}{$rep}[$lane-1] = $count;
  }
  return %h
}

sub GetRep
{
  my ($cond, $rep) = @_;
  
  my $file = CountFile(undef, $cond, $rep, $type);
  my ($d) = ReadExpressionFile($file, nonzero=>0);
  #my $di = floor($d + 0.5);
  
  my $N = sumover($d>0);
  my $sel = which($N >= 1);
  $d = $d(,$sel);
  
  unless($type eq 'lev')
  {
    my $tot = sumover(transpose($d));   # sum of mapped reads
    $tot /= sum($tot) / nelem($tot);
  
    my $ltot = pdl(@{$lanetot{$cond}{$rep}});  # sum of all reads
    $ltot /= sum($ltot) / nelem($ltot);
  
    #print "$cond $rep $tot $ltot\n";
  
    $d /= $ltot;                           # normalize
    #$d = floor($d + 0.5);                 # need integers
  }
  
  #my ($ncol, $nrow) = dims($d);
  #print "$ncol x $nrow\n";
  
  return $d
}



sub PlotPDist
{
  my ($w, $pan, $P, $ns, $limit) = @_;
  
  my $z = 1e-9;
  $P->where($P <= $z) .= $z;
      
  my ($x, $h) = BuildHistogram(log10($P), undef, undef, undef, 1, 1);
  
  PlotPanelN($w, $pan,
    BOX => [log10($z)-0.5, 0, 0, 1.99],
    XLAB => 'log P',
    YLAB => 'Frequency',
    CHARSIZE => 0.3,
    MAJTICKSIZE=>0.2,
    MINTICKSIZE=>0.2,
    xbox => 'I', ybox => 'I'
  );
  
  BarPlot($w, $x, $h, colour=>[215,215,255]);
  HistPlot($w, $x, $h);
  
  my $l = log10($limit);
  LinePlot($w, pdl($l,$l), pdl(0,10), COLOR=>'RED');
  
  TopText($w, $pan, x=>0.1, CHARSIZE=>0.8);
  
  #my $s = sprintf "%3.1f %.2g", 100*$n1/$n, 100*$n2/$n;
  my $s = sprintf "n#ds#u=%d", $ns;
  TopText($w, $s, x=>0.1, y=>-5, CHARSIZE=>0.4);
}




=head1 SYNOPSIS

  poisson_reps.pl  
    
=head1 OPTIONS

=over 5

=item B<-readcountfile>=I<pathname>

The file with total read counts per lane, created with C<count_fastq_reads_lane.pl>. The default name is C<all_readcount_lane.dat>.

=item B<-outfile>=I<pathname>

The file with test results. The default name is C<all_poiss.dat>. There are four columns written to this file: mean, standard deviation, p-value (not adjusted) and (chi-square distributed) T-statistic.

=item B<-ptype>=[dispersion|logratio]

Type of the test.


=back
