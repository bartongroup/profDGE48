#!/sw/bin/perl

=head1 NAME

B<compare_null.pl>

=head1 DESCRIPTION

Compare results from null test - FPR bootstraps for the same condition. Requires a test description file with all details of the bootstraps (de_tests_samecond.txt is default).

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;
use Switch 'Perl6';

use GRNASeq;
use GetOpt;
use Powertests;

use PLgraphs;
use Tools;
use Stats;
use Distribution;

$| = 1;
Stamp();

#my $multicor = 'bh';

my $graph = 'all';
my $testinfofile = "de_tests_samecond.txt";
my $dir = 'powerstats_same_db';
my $tests = 'lt,bayseq,cuffdiff,deseq1,ebseq,edger,edgerglm,limma,poissonseq,samseq';
my $maxy = 0.4;
my ($help, $man);
GetMyOptions(
  'dir=s' => \$dir,
  'maxy=f' => \$maxy,
  'testinfofile=s' => \$testinfofile,
  'tests=s' => \$tests,
  'graph=s' => \$graph,
   help => \$help,
   man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;

my $cs = 0.45;
my $cs2 = 0.3;
my @nrep = (3, 8, 20);
my @col = @contrast_colours;

###############################################

my %tests = ReadTestfileInfo($testinfofile);
my @tests = split(/\,/, $tests);
my $N = @tests;

my $w = NewWindow($psfile);

given($graph)
{
  when 'pool' {CompareNsigPool()}
  when 'prop' {PlotNplus()}
  when 'all' {PlotNsig()}
}

$w->close();


##########################################

sub PlotNplus
{
  my $div = 0.95;
  
  $PL_ymin = 0.08;
  $PL_xmin = 0.08;
  $PL_xmax = $div;
  
  $PL_nx = 4;
  $PL_ny = 3;
  $PL_deltax = 0.03;
  $PL_deltay = 0.06;
  
  
  my $i = 0;
  for my $test (@tests)
  {
    PlotPanelN($w, $i+1,
      BOX => [0, 21, 0, 1],
      XLAB => '## replicates',
      #YLAB => 'Bootstraps with FP fraction > 0.05'
      YLAB => 'Prop(FPR > 0.05)',
      xbox => 'I', ybox => 'I'
    );
    
    #my $col = $contrast_colours[$i];
    my $col = 'BLACK';
    
    TopText($w, $tests{$tests[$i]}{name}, COLOR=>$col, y=>1);
    
    my @nr = ();
    my @frac = ();
    my @prop = ();
    my @err = ();
    for my $nrep (3 .. 20)
    {
      my $d = GetNsig($test, $nrep);
      next unless defined $d;
      my $s = nelem($d->where($d>0.05));
      my $n = nelem($d);
      my ($frac, $err) = PropCI($n, $s);
      push @nr, $nrep;
      push @frac, $frac;
      push @prop, $s/$n;
      push @err, $err;
    }
    my $nr = pdl(@nr);
    my $frac = pdl(@frac);
    my $err = pdl(@err);
    my $prop = pdl(@prop);
    
    #wcols $nr, $frac, $err;
    
    #LinePlot($w, $nr, $frac, COLOR=>$contrast_colours[$i], LINEWIDTH=>3);
    #LinePlot($w, $nr, $frac, COLOR=>$col, LINEWIDTH=>1);
    PointPlot($w, $nr, $frac, COLOR=>$col, SYMBOLSIZE=>0.001, YERRORBAR=>2*$err);
    PointPlot($w, $nr, $prop, COLOR=>$col, SYMBOLSIZE=>2);
    #HistPlot($w, $nr, $frac, colour=>$contrast_colours[$i], LINEWIDTH=>3);
    $i++; 
  }
  
  #Legend($div);
  
}

##########################################

sub PlotNsig
{
  my $div = 0.95;
  
  $PL_ymin = 0.08;
  $PL_xmin = 0.08;
  $PL_xmax = $div;
  
  $PL_nx = 4;
  $PL_ny = 3;
  $PL_deltax = 0.03;
  $PL_deltay = 0.06;
  
  
  my $i = 0;
  for my $test (@tests)
  {
    my $mm = ($test eq 'degseq') ? 1 : $maxy;
    my $ybox = ($test eq 'degseq') ? 'IN' : 'I';
    my $pan = ($test eq 'bobfish') ? 12 : $i + 1;
    PlotPanelN($w, $pan,
      BOX => [0, 21, 0, $mm],
      XLAB => '## replicates',
      #YLAB => 'Bootstraps with FP fraction > 0.05'
      YLAB => 'FPR',
      xbox => 'I', ybox => $ybox
    );
    
    #my $col = $contrast_colours[$i];
    my $col = 'BLACK';
    
    TopText($w, $tests{$tests[$i]}{name}, COLOR=>$col, y=>1);
    
    my @nr = ();
    my @frac = ();
    my @prop = ();
    my @err = ();
    my @d = ();
    for my $nrep (3 .. 20)
    {
      my $d = GetNsig($test, $nrep, $dir);
      next unless defined $d;
      my $s = nelem($d->where($d>0.05));
      my $n = nelem($d);
      my ($frac, $err) = PropCI($n, $s);
      push @nr, $nrep;
      push @frac, $frac;
      push @prop, $s/$n;
      push @err, $err;
      push @d, $d;
    }
    my $nr = pdl(@nr);
    my $frac = pdl(@frac);
    my $err = pdl(@err);
    my $prop = pdl(@prop);
    
    #wcols $nr, $frac, $err;
    
    #LinePlot($w, $nr, $frac, COLOR=>$contrast_colours[$i], LINEWIDTH=>3);
    #LinePlot($w, $nr, $frac, COLOR=>$col, LINEWIDTH=>1);
    #PointPlot($w, $nr, $frac, COLOR=>$col, SYMBOLSIZE=>0.001, YERRORBAR=>2*$err);
    #PointPlot($w, $nr, $prop, COLOR=>$col, SYMBOLSIZE=>2);
    #HistPlot($w, $nr, $frac, colour=>$contrast_colours[$i], LINEWIDTH=>3);
    BoxPlot($w, \@d, solidwhisk=>1, colour=>[215,215,255], outliers=>1, pos=>3);
    LinePlot($w, pdl(0,21), pdl(0.05,0.05), COLOR=>'RED');
    
    $i++; 
  }
  
  #Legend($div);
  
}

##########################################

sub CompareNsigPool
{
  
  my $div = 0.2;
  
  $PL_ymin = 0.6;
  $PL_ymax = 0.9;
  $PL_xmin = 0.05;
  $PL_xmax = $div;
  PlotPanelN($w, 1, %PL_empty, BOX => [0,1,0,$N+1]);
  for my $i (1 .. $N)
  {
    my $name = $tests{$tests[$i-1]}{name};
    $w->text($name, COLOR=>'BLACK', CHARSIZE=>0.7, TEXTPOSITION=>[0.9, $N-$i+1, 0, 0, 1])
  }
  
  
  $PL_nx = 1;
  $PL_ny = 1;
  $PL_deltax = 0.06;
  $PL_xmin = $div;
  $PL_xmax = 0.70;
  
  $maxy = 1;
  PlotNsigNull(1, 0);
  
  #FullBox($w);
  #TopText($w, "Null test", c=>0.5, x=>0.35, y=>-2);
}


##########################################

sub CompareNsig
{
  
  my $div = 0.2;
  
  $PL_ymin = 0.6;
  $PL_ymax = 0.9;
  $PL_ymax = 0.9;
  $PL_xmin = 0.05;
  $PL_xmax = $div;
  PlotPanelN($w, 1, %PL_empty, BOX => [0,1,0,$N+1]);
  for my $i (1 .. $N)
  {
    my $name = $tests{$tests[$i-1]}{name};
    $w->text($name, COLOR=>'BLACK', CHARSIZE=>0.7, TEXTPOSITION=>[0.9, $i, 0, 0, 1])
  }
  
  
  $PL_nx = @nrep;
  $PL_ny = 1;
  $PL_deltax = 0.06;
  $PL_xmin = $div;
  $PL_xmax = 0.95;
  
  
  my $pan = 1;
  
  for my $nrep (@nrep)
    {PlotNsigNull($pan++, $nrep)}
  
  FullBox($w);
  TopText($w, "Null test", c=>0.5, x=>0.35, y=>-2);
}

##########################################

sub PlotNsigNull
{
  my ($pan, $nrep) = @_;  
    
  PlotPanelN($w, $pan,
      BOX => [0, $maxy, -1, $N],
      XBOX => 'BINST',
      YBOX => '',
      XLAB => 'FP fraction',
      #XTICK => 0.1, NXSUB=>2
    );
    
  my @dat = ();
  local *F;
  
  for my $test (@tests)
  {
    my $d;
    if($nrep == 0)
    {
      my @d = ();
      for my $r (3 .. 20)
      {
        my $dd = GetNsig($test, $r);
        #print nelem($dd), " ";
        push @d, $dd;
      }
      $d = pdl(@d)->flat();
      print nelem($d), "\n";
    }
    else
      {$d = GetNsig($test, $nrep)}
    push @dat, $d;
  }
  
  @dat = reverse @dat;
  
  BoxPlot($w, \@dat, swap=>1, pcolour=>[110,110,110], outliers=>1, solidwhisk=>1, colour=>[215,215,255],);
  LinePlot($w, pdl(0.05,0.05), pdl(-1,$N), COLOR=>'RED');
}

###########################################

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

###########################################

sub GetNsig
{
  my ($test, $nrep, $dir) = @_;
  
  local *F;
  my $file = "$dir/power_${test}_true-${test}_nsig.stats";
  return unless -e $file;
  open F, $file or die "Cannot open $file\n";
  <F>;
  my $d;
  while(my $line = <F>)
  {
    my @d = split /\t/, $line;
    my $n = shift @d;
    next unless $n == $nrep;
    $d = pdl(@d);
  }
  
  return $d;
}

###########################################

sub Legend
{
  my ($div) = @_;
  
  $PL_xmin = $div;
  $PL_xmax = 0.99;
  $PL_nx = 1;
  $PL_ny = 1;
  PlotPanelN($w, 1, %PL_empty);
  
  my $i = 0;
  for my $i (0 .. @tests - 1)
  {
    my $x = 0.1;
    my $y = 0.9 - $i * 0.04;
    my $dx = 0.2;
    
    my $c = $col[$i];
 
    LinePlot($w, pdl($x, $x+$dx), pdl($y, $y), COLOR=>$c, LINEWIDTH=>3);
    $w->text($tests{$tests[$i]}{name}, COLOR=>$c, CHARSIZE=>0.5, TEXTPOSITION=>[$x+$dx+0.05,  $y, 0, 0, 0]);
  }
}  



=head1 SYNOPSIS

  compare_null.ps -psfile=null.ps

=head1 OPTIONS

=over 4

=item B<-testinfofile>=I<pathname>

A file with DE tools metadata (for the 'same condition' bootstraps). See C<plot_fcp.pl> for details. Default name is C<de_tests_samecond.txt>

=item B<-dir>=I<string>

Directory where the results from C<make_powerstats_db.pl> are stored. The default value is C<powerstats_same_db>.

=item B<-tests>=I<string>

A comma-delimited list of tests/tools to be included in the plot. The default value is 'lt,bayseq,cuffdiff,degseq,deseq1,deseq2,ebseq,edger,limma,noiseq,samseq,poissonseq'.

=item B<-graph>=I<string>

Type of graph to plot. The default value is 'all'. Do not change!

=back
