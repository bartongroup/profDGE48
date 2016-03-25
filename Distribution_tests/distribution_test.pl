#!/sw/bin/perl

=head1 NAME

B<distribution_test.pl>

=head1 DESCRIPTION

Performs goodness-of-fit distribution test for gene expression data across biological replicates. For normal, log-normal and Poisson distributions tests are performed by this script. For across-lane Poisson and negative binomial, there are separate scripts, and this scritp only reads data.

A figure is created with test results, showing distribution of p-values and p-value versus mean count plot.

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

use CompBio::Normality;
use CompBio::Distribution;
use CompBio::Statistics;
use CompBio::PLgraphs;
use CompBio::Tools;
use CompBio::Heatmap;
use Switch 'Perl6';

Stamp();
$| = 1;

my %distribution = (
  norm => 'normal',
  lnorm => 'lognormal',
  pois => 'Poisson',
  nb => 'negative binomial (M)',
  nbks => 'negative binomial (KS)',
  nbhg => 'negative binomial (HG)'
);

my %testtype = (
  dispersion => 'dispersion',
  logratio => 'log ratio'
);

ReadDefs();

my $lpoisfile = 'all_poiss.dat';
my $dist = 'pois';
my $ptype = 'dispersion';
my $limit = 0.05;
my $out;
my $sigma = 5;
my $nbstat = 'm';
my ($withsim, $myminp, $report);
my ($help, $man);
GetMyOptions(
  'lpoisfile=s' => \$lpoisfile,
  'dist=s' => \$dist,
  'ptype=s' => \$ptype,
  'minp=f' => \$myminp,
  'out=s' => \$out,
  'sigma=f' => \$sigma,
  'nbstat=s' => \$nbstat,
  withsim => \$withsim,
  report => \$report,
);

$nonzero = -1 if $dist eq 'lnorm';

my $zero = 1e-17;

my ($minm, $maxm) = (-2, 6);
my ($minvm, $maxvm) = (-1, 5);
my ($minp, $maxp) = (-17.1, 0);
$minp = $myminp if defined $myminp; 


my ($d, $di, $m, $s, $nrep, $ngen, $genes, $name, $sel, $P);

if($dist eq 'nb')
{
  my $cl = ($clean) ? '_clean' : '';
  $cl = '_spclean' if $spclean;
  my $tp = ($type eq 'lev') ? "_lev" : '';
  my $file = "${cond}_nbtest_${nbstat}${tp}${cl}_${norm}.dat";
  #my $file = "${cond}_nbtest${tp}${cl}.dat";
  print "Reading $file...";
  ($m, $P) = rcols $file, 1, 3;
  print " done\n";
  $ngen = nelem($m);
  $name = "${cond}${tp}${cl}";
  $nonzero = 1;
  
  my @g = ReadGeneList($file);
  $genes = \@g;
  $sel = sequence(scalar @g);
  
  #undef $report;
}
elsif($dist eq 'tpois')
{
  die "Cannot find file $lpoisfile\n" unless -e $lpoisfile;
  ($m, $P) = rcols $lpoisfile, 0, 2;
  $ngen = nelem($m);
  $name = 'Technical Poisson';
}
else
{
  ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel) = GetNormalizedData(
    condition=>$cond,
    replicate=>$rep,
    exclude=>$exclude,
    nonzero=>$nonzero,
    clean=>$clean,
    spclean=>$spclean,
    norm=>$norm,
    clip=>$clip,
    type=>$type
  );
  $di = floor($d + 0.5);
  print "$ngen \(nz = $nonzero) genes in $nrep replicates read.\n";
}

my $sim;
if($withsim)
{
  print "Generating simulated data...       ";
  $sim = SimulatedData();
  print "\b\b\b\b\b\b\b done  \n";
}

print "Testing\n";
my $Ps;
if($dist eq 'norm')
{
  $P = NormalityFromData($d, 0, $out, $sigma);
  $Ps = NormalityFromData($sim) if $withsim;
}
elsif($dist eq 'lnorm')
{
  $P = NormalityFromData($d, 1, $out, $sigma);
  $Ps = NormalityFromData($sim, 1) if $withsim;
}
elsif($dist eq 'pois')
{
  $P = PoissonFromData($di, $ptype);
  $Ps = PoissonFromData($sim, $ptype) if $withsim;
}
elsif($dist eq 'nbks')
{
  $P = NBKSFromData($di);
}
elsif($dist eq 'nbhg')
{
  my $N = sumover($di==0);
  my $zel = which($N > 0);
  
  #$sel = which($m < 1e9);
  
  my $d0 = $di(,$zel);
  $m = $m($zel);
  $ngen = nelem($zel);
  print "$ngen with zeroes\n";
  $P = NBHGFromData($d0);
  
  $sel = $sel($zel);
  
  #$Ps = NBFromData($sim, $ptype) if $withsim;
}



print "done\n";

#wcols $m, $s, $s**2/$m, $P;
my $cl = ($clean) ? '_clean' : '';
$cl = '_spclean' if $spclean;
#my $outfile = "${cond}_${dist}_${norm}${cl}_${multicor}_test.dat";
my $outfile = "${cond}_${dist}_${norm}${cl}_test.dat";
ReportGenes($outfile, $genes, $sel, $m, $P) if $report;
  
my $win = NewWindow($psfile);
$PL_nx = 2;
$PL_ny = 2;
$PL_deltay = 0.06;
$PL_deltax = 0.06;

my $pan = 1;

PlotPDist($pan++, $P);
#PlotPMean($pan++, $m, $P, '');
PlotSigfracMean($pan++, $m, $P);

#PlotP(2, $sim, $Ps, 'Simulated') if $withsim;

#PlotPDist(2, $P);
#PlotPDist(2, $sP);

my $str = sprintf "Test for %s distribution (nz = %d, n = %d) norm=%s mult=%s", $distribution{$dist}, $nonzero, $ngen, $norm, $multicor;
$str .= sprintf " (%s test)", $testtype{$ptype} if $dist eq 'pois';
FullBox($win);
TopText($win, $str, x=>0.07, y=>-1, CHARSIZE => 0.65);
TopText($win, $name, x=>0.07, y=>-3, CHARSIZE => 0.65);

$win->close();

########################################################

sub PlotPMean
{
  my ($pan, $m, $P, $lab) = @_;


  my ($Plim, $ns);
  if($multicor eq 'hb')
   {($Plim, $ns) = HolmBonferroniLimit($P, $limit)}
  elsif($multicor eq 'bh')
    {($Plim, $ns) = BenjaminiHochbergLimit($P, $limit)}
  else
    {die "Huh?\n"}

  $P->where($P==0) .= $zero;
  my $y = log10($P);
  my $x = log10($m);

  $lab ||= '';
  my ($minx, $maxx) = BoxMinMax($x);
  my ($miny, $maxy) = BoxMinMax($y);
  #print "$maxy\n";
  
  PlotPanelN($win, $pan,
    BOX => [$minm, $maxm, $minp, $maxp],
    XLAB => 'log mean',
    YLAB => 'log p',
    xbox=>'NI', ybox=>'NI',
    forceylab=>1,
    forcexlab=>1,
  );
  
  my $sel1 = which($P>$zero);
  my $sel2 = which($P<=$zero);
  PointPlot($win, $x($sel1), $y($sel1));
  PointPlot($win, $x($sel2), $y($sel2), COLOR=>'GOLD2');

  TopText($win, $lab, y=>1);
  
  #limline($y, $limit, $ngen);
  #limline($y, $limit/$ngen, $ngen);
  limline($y, $Plim, $ngen) if $Plim > 0;
}

########################################################

sub PlotSigfracMean
{
  my ($pan, $m, $P, $lab) = @_;


  my ($Plim, $ns);
  if($multicor eq 'hb')
   {($Plim, $ns) = HolmBonferroniLimit($P, $limit)}
  elsif($multicor eq 'bh')
    {($Plim, $ns) = BenjaminiHochbergLimit($P, $limit)}
  else
    {die "Huh?\n"}

  $P->where($P==0) .= $zero;
  my $f = $P <= 0.05;
  my $lm = log10($m);
  
  my $nbin = 30;
  my ($min, $max) = (-1, 7);
  my $step = ($max - $min) / 30;
  my @x = ();
  my @y = ();
  for my $i (1 .. $nbin)
  {
    my $x = $min + $step * ($i - 1);
    my $x1 = $x - $step/2;
    my $x2 = $x + $step/2;
    my $ff = $f->where(($lm >= $x1) & ($lm < $x2));
    my $frac = 0;
    $frac = sum($ff) / nelem($ff) if nelem($ff) > 0;
    push @x, $x;
    push @y, $frac
  } 

 
  PlotPanelN($win, $pan,
    BOX => [$min, $max, 0, 1],
    XLAB => 'log mean',
    YLAB => 'Significan fraction',
    xbox=>'NI', ybox=>'NI',
    forceylab=>1,
    forcexlab=>1,
  );
  
  HistPlot($win, pdl(@x), pdl(@y));

}

sub limline
{
  my ($P, $lim, $n) = @_;
  
  my $q = log10($lim);
  my $out = nelem(which($P < $q));
  #print "out = ", 100 * $out / $n, "\n";
  my $str = sprintf "%.2g%%",100 * $out / $n;
  LinePlot($win, pdl($minm,$maxm),  pdl($q,$q), COLOR=>'RED');
  my $xx = $maxm - 0.05*($maxm-$minm);
  #$win->text($str, COLOR=>'RED', CHARSIZE=>0.6, TEXTPOSITION=>[$xx, $q+0.15, 0, 0, 0]);  
}

########################################################

sub ReportGenes
{
  my ($file, $genes, $sel, @val) = @_;
  
  local *F;
  open F, ">$file" or die;
  for my $i (0 .. nelem($sel) - 1)
  {
    my $gene = $genes->[at($sel, $i)];
    my @s = ($gene);
    for my $v (@val)
      {push @s, at($v, $i)}
    my $line = join "\t", @s;
    print F "$line\n";
  }
  print "File $file created.\n"
}

########################################################

sub PlotPDist
{
  my ($pan, $P) = @_;

  my ($Plim, $ns);
  if($multicor eq 'hb')
   {($Plim, $ns) = HolmBonferroniLimit($P, $limit)}
  elsif($multicor eq 'bh')
    {($Plim, $ns) = BenjaminiHochbergLimit($P, $limit)}
  else
    {die "Huh?\n"}
  
  print "Plim = $Plim\n";
  $Plim = ($Plim > 0) ? log10($Plim) : -1;

  my ($x, $h) = BuildHistogram(log10($P->where($P)), undef, undef, undef, 1, 1);
  my $max = max($h)*1.03;

  PlotPanelN($win, $pan,
    BOX => [-17, 0, 0, $max],
    %PL_empty
  );

  BarPlot($win, $x, $h, colour=>[215,215,255]);
  HistPlot($win, $x, $h);
  LinePlot($win, pdl($Plim, $Plim), pdl(0, $max), COLOR=>'RED');

  PlotPanelN($win, $pan,
    BOX => [-17, 0, 0, max($h)*1.03],
    XLAB => 'log p',
    YLAB => 'Frequency',
    xbox => 'NI', ybox => 'I',
    forceylab=>1,
    forcexlab=>1,
  );
  
  my $s = sprintf "n#ds#u=%d (%5.2f%%)", $ns, 100*$ns/$ngen;
  TopText($win, $s, x=>0.05, y=>-1, CHARSIZE=>0.8);
}

sub PlotSMDist
{
  my ($pan, $r) = @_;
  
  my $lr = log10($r);
  my ($x, $h) = BuildHistogram($lr, undef, undef, undef, 1);
  
  PlotPanelN($win, $pan,
    BOX => [min($x), max($x), 0, max($h)*1.03],
    XLAB => 'log S#u2#d/M',
    YLAB => 'Normalized frequency',
  );
  
  BarPlot($win, $x, $h, colour=>[215, 215, 255]);
  HistPlot($win, $x, $h, linewidth=>4);
  OverPlot($win, pdl(0,0), pdl(0,10), PLOTTYPE=>'LINE', COLOR=>'RED', LINESTYLE=>4, LINEWIDTH=>2);}



sub SimulatedData
{
  my @dat = ();
  for my $i (0 .. $ngen-1)
  { 
    Percent1($i/$ngen) if $i % 100 == 0;
    my $M = at($m, $i);
    my $S = at($s, $i);
    
    my $x;
    if($dist eq 'norm')
      {$x = GenerateLognormal($M, $S, $nrep)}
    elsif($dist eq 'pois')
      {$x = GeneratePoisson($M, $nrep)}
    push @dat, $x;
  }
  my $dat = pdl(@dat);
  
  return $dat;  
}



=head1 SYNOPSIS

  distribution_test.pl -dist=lnorm -clean -psfile=lnrom_test.ps
    
=head1 OPTIONS

=over 5

=item B<-dist>=I<string>

Distribution to use for the test.

  norm - normal
  lnorm - log-normal
  pois - Poisson
  tpois - technical Poisson (between lanes); data are read from "all_poiss.dat" file, created with poisson_reps.pl
  nb - negative binomial using Meintanis method; you need to run grid_launcher_nbtest.pl first

Other NB methods do not work very well.

=item B<-lpoisfile>=<pathname>

File with lane goodness-of-fit results, created by C<poisson_reps.pl>.

=item B<-clean>

If defined, only clean replicates will be used, as defined in C<defs.dat> file.

=item B<-norm>=I<[none|totcount|totcountm|deseq|tmm|spikein]>

Normalization of gene expression data.

  * none - raw counts used (not recommended)
  * totcount - total count divided by 1e6 (not recommended)
  * totcountm - total count devided by its mean
  * deseq - DESeq normalization
  * tmm - TMM normalization
  * spikein - normalization to spike-ins (requires files with spike-in normalizing factors, <condition>_spikein_normfac.txt in countsdir directory, as defined in defs.dat)

=item B<-ptype>=I<[dispersion|logratio]>

Type of Poisson test used. Default is C<dispersion>.

=item B<-minp>=I<number>

Minimum p-value in the plot.

=item B<-out>=I<[chauvenet|sigma]>

The method to reject outliers before doing test for normality or log-normality. If not specified, outliers are not rejected.

=item B<-sigma>=I<number>

Sigma parameter for C<sigma> method of rejecting outliers.

=item B<-nbstat>=I<[m|a]>

Statistic to use with negative binomial test. 'm' - Meintanis, 'a' - Anderson-Darling. The latter one doesn't work very well. Default value is 'm'.

=item B<-report>

If specified a file will be created with mean and p-values, which can used for further analysis. 

=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut