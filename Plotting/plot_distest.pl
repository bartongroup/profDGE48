#!/sw/bin/perl

=head1 NAME

B<plot_distest.pl>

=head1 DESCRIPTION

Plot a figure with results from log-normal and negative binomial tests. Need to run C<distribution_test.pl> and C<grid_launcher_nb_test.pl> first to create test results. These test results should be in files:

  WT_lnorm_lev_clean_test.dat
  Snf2_lnorm_lev_clean_test.dat
  WT_lnorm_lev_test.dat
  Snf2_lnorm_lev_test.dat
  WT_nb_deseq_lev_test.dat
  Snf2_nb_deseq_lev_test.dat
  WT_nb_lev_test.dat
  Snf2_nb_lev_test.dat

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use GetOpt;

use Distribution;
use Stats;
use PLgraphs;
use Tools;

Stamp();

#my $cond = 'WT';
my $limit = 0.05;
my $zero = 0.2e-17;
my $psfile;
my $norm = 'lev';
my $multicor = 'bh';
my ($help, $man);
GetOptions(
 # 'cond=s' => \$cond,
 'norm=s' => \$norm,
  'psfile=s' => \$psfile,
  #'multicor=s' => \$multicor,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;




my $w = NewWindow($psfile);

$PL_nx = 4;
$PL_ny = 3;
$PL_deltax = 0.02;
$PL_deltay = 0.05;
$PL_ymin = 0.2;

my ($minm, $maxm) = (-1.9999, 6);
my ($minp, $maxp) = (-18, 0);


my $pan = 1;
for my $test qw(norm lnorm nb)
{
  for my $cl ('_clean', '')
  {
    for my $cond qw(WT Snf2)
    { 
      #my $file = "${cond}_${test}_${norm}${cl}_${multicor}_test.dat";
      my $file = "${cond}_${test}_${norm}${cl}_test.dat";
      next unless -e $file;
      print "$file\n";
      
      my $ppos = ($test eq 'nb' && $norm eq 'lev') ? 3 : 2;  # old format
      my ($m, $P) = rcols $file, 1, $ppos;
      PlotPMean($w, $pan, $m, $P);
      my $c = ($cl eq '_clean') ? 'clean' : '';
      my $t;
      if($test eq 'norm') {$t = 'nm'}
      elsif($test eq 'lnorm') {$t = 'ln'}
      else {$t = 'nb'}

      my $l = chr($pan+96);
      TopText($w, "\($l) $t $cond $c", CHARSIZE=>0.7, x=>0, y=>1, COLOR=>[110,110,110]);
      $pan++
    }
  }
}

$w->close();


sub PlotPMean
{
  my ($win, $pan, $m, $P, $lab) = @_;

  my $ngen = nelem($m);

  my ($Plim, $ns) = HolmBonferroniLimit($P, $limit);

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
    xbox=>'I', ybox=>'I',
    #forceylab=>1,
    #forcexlab=>1,
  );
  
  my $sel1 = which($P>$zero);
  my $sel2 = which($P<=$zero);
  PointPlot($win, $x($sel1), $y($sel1));
  PointPlot($win, $x($sel2), $y($sel2), COLOR=>'GOLD2');

  #TopText($win, $lab, y=>1);
  $win->text("$ns/$ngen", CHARSIZE=>0.4, COLOR=>'BLUE', TEXTPOSITION=>[$maxm-0.0, -16.8, 0, 0, 1]);
  
  #limline($y, $limit, $ngen);
  #limline($y, $limit/$ngen, $ngen);
  limline($win, $y, $Plim, $ngen) if $Plim > 0;
}

sub limline
{
  my ($win, $P, $lim, $n) = @_;
  
  my $q = log10($lim);
  my $out = nelem(which($P < $q));
  #print "out = ", 100 * $out / $n, "\n";
  my $str = sprintf "%.2g%%",100 * $out / $n;
  LinePlot($win, pdl($minm,$maxm),  pdl($q,$q), COLOR=>'RED');
  my $xx = $maxm - 0.05*($maxm-$minm);
  #$win->text($str, COLOR=>'RED', CHARSIZE=>0.6, TEXTPOSITION=>[$xx, $q+0.15, 0, 0, 0]);  
}




=head1 SYNOPSIS

  plot_distest.pl -norm=lev -psfile=dist_tests.ps
    
=head1 OPTIONS

=item B<-norm>=I<string>

Normalization used in distribution test scripts (see DESCRIPTION). This is used to identify the test result files. The default value is 'lev'.

=item B<-psfile>=I<pathname>

Output postscript file.

=cut
