package PLgraphs;
require Exporter;
require PDL;
require PDL::NiceSlice;
require PDL::Graphics::PLplot;
require PDL::Options;
require Heatmap;
require Stats;
require Carp;

=head1 NAME

CompBio::PLgraphs - graphical routines based on PLplot

=head1 SYNOPSIS

  use CompBio::PLgraph;
  my ($boxmin, $boxmax) = BoxMinMax($x, $spacer);
  my $win = NewWindow();
  BoxPlot($win, $data, colour => [100, 200, 200], linewidth => 2, cloud => 1);
  my ($x1, $x2, $h) = BuildHistogram2($d, $n, $min, $max, $renorm);
  my ($x, $h) = BuildHistogram($d, $n, $min, $max, $renorm);
  BarPlot($win, $x, $y, colour => [100, 200, 200], width => 0.8, boxes => 1)
  BarPlot2($win, $x1, $y, y0 => 0);
  HistPlot($win, $x1, $x2, $y, colour => 'RED', linewidth => 2, linestyle => 1);
  my $colour = GraphColour($i);

=head1 FUNCTIONS

=over 4

=cut

use strict;
use warnings;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;
use PDL::Options;
use Heatmap;
use Stats qw(corr TrimmedMean);
use Carp;

our $VERSION = 2.2;

our @ISA = ("Exporter");
our @EXPORT = qw(
  NewWindow
  PanelCoord
  PanelCoordN
  PlotPanelXY
  PlotPanelN
  OverPlot
  LinePlot
  ShadeErrPlot
  XYShadePlot
  MAPlot
  PointPlot
  TopText
  FullBox
  KDE
  BuildHistogram
  BuildHistogram2
  MultiBarPlot
  BarPlot
  BarPlot2
  HistPlot
  HistPlot2
  BoxPlot
  DistPlot
  BoxMinMax
  GraphColour
  percentile
  PLStore
  PLRecover
  $PL_nx
  $PL_ny
  $PL_xmin
  $PL_xmax
  $PL_ymin
  $PL_ymax
  $PL_deltax
  $PL_deltay
  %PL_empty
  @contrast_colours
);  

our $PL_nx = 1;
our $PL_ny = 1;
our $PL_xmin = 0.08;
our $PL_xmax = 0.95;
our $PL_ymin = 0.08;
our $PL_ymax = 0.95;
our $PL_deltax = 0.01;
our $PL_deltay = 0.01;
our %PL_empty = (XBOX => '', YBOX => '', XLAB => '', YLAB => '');

sub cc2t { [map {hex} split ' ', shift] }
my %colours = (
               BLACK          => [  0,  0,  0],
               RED            => [240, 50, 50],
               YELLOW         => [255,255,  0],
               GREEN          => [  0,255,  0],
               AQUAMARINE     => [127,255,212],
               PINK           => [255,192,203],
               WHEAT          => [245,222,179],
               GREY           => [190,190,190],
               BROWN          => [165, 42, 42],
               BLUE           => [  0,  0,255],
               BLUEVIOLET     => [138, 43,226],
               CYAN           => [  0,255,255],
               TURQUOISE      => [ 64,224,208],
               MAGENTA        => [255,  0,255],
               SALMON         => [250,128,114],
               WHITE          => [255,255,255],
               ROYALBLUE      => cc2t('2B 60 DE'),
               DEEPSKYBLUE    => cc2t('3B B9 FF'),
               VIOLET         => cc2t('8D 38 C9'),
               STEELBLUE1     => cc2t('5C B3 FF'),
               DEEPPINK       => cc2t('F5 28 87'),
               MAGENTA        => cc2t('FF 00 FF'),
               DARKORCHID1    => cc2t('B0 41 FF'),
               PALEVIOLETRED2 => cc2t('E5 6E 94'),
               TURQUOISE1     => cc2t('52 F3 FF'),
               LIGHTSEAGREEN  => cc2t('3E A9 9F'),
               SKYBLUE        => cc2t('66 98 FF'),
               FORESTGREEN    => cc2t('4E 92 58'),
               CHARTREUSE3    => cc2t('6C C4 17'),
               GOLD2          => cc2t('EA C1 17'),
               SIENNA1        => cc2t('F8 74 31'),
               CORAL          => cc2t('F7 65 41'),
               HOTPINK        => cc2t('F6 60 AB'),
               LIGHTCORAL     => cc2t('E7 74 71'),
               LIGHTPINK1     => cc2t('F9 A7 B0'),
               LIGHTGOLDENROD => cc2t('EC D8 72'),
);

our @contrast_colours = (
  [0, 0, 255],
  [0, 127, 0],
  [255, 0, 0],
  [0, 191, 191],
  [191, 0, 191],
  [191, 191, 0],
  [63, 63, 63],
  [191, 63, 63],
  [242, 242, 0],
  [63, 63, 191],
  [191, 191, 191],
  [0, 255, 0],
  [193, 145, 43],
  [137, 160, 56],
  [86, 145, 234],
  [255, 25, 153],
  [224, 191, 186],
  [25, 124, 119],
  [168, 86, 165],
  [252, 104, 58]
);

my $PI = 4 * atan2(1, 1);

#########################################################

=item NewWindow

  my $win = NewWindow();
  my $win = NewWindow($psfile);
  
Opens and returns PLplot window. If $psfile is not specified, opens an X-windows.

=cut

sub NewWindow
{
  my ($file, $orientation) = @_;
  
  my $inch = 2.54;
  my $DPI = 64;
  my $A4_x = 21.0 / $inch * $DPI;
  my $A4_y = 29.7 / $inch * $DPI; 
  my $win;
  
  if(defined $file)
  {
    if($file =~ /\.ps$/)
    {
      $orientation = 1 unless defined $orientation;
      $win = PDL::Graphics::PLplot->new(
        #DEV => 'psc',
        DEV => 'psc',
        FILE => $file,
        PAGESIZE => [$A4_x, $A4_y],
        ORIENTATION => $orientation,
        #OPTS => {drvopt=>'hrshsym=1'}
      )
    }
    elsif($file =~ /\.svg$/)
    {
      $win = PDL::Graphics::PLplot->new(
        DEV => 'svg',
        FILE => $file,
        PAGESIZE => [1024, 1024],
      )
    }
  }
  else
  {
    $orientation = 0 unless defined $orientation;
    $win = PDL::Graphics::PLplot->new(
      DEV => 'xwin',
      PAGESIZE => [900, 900],
      ORIENTATION => $orientation,
    );
  }
  return $win;
}

####################################################

=item PanelCoord

  my ($x1, $x2, $y1, $y2) = PanelCoord($x, $y);

Returns coordinates (to be used in PlotPosition inside env) of a panel ($x, $y). There are some global variables that can be changed:

  $PL_nx, $PL_ny - number of panels in each direction
  $PL_xmin, $PL_xmax, $PL_ymin, $PL_ymax - borders of the entire plot (between 0 and 1)
  $PL_deltax, $PL_deltay - gaps between panels (in 0-1 units)

=cut

sub PanelCoord
{
  my ($x, $y) = @_;
  
  my $sizex = ($PL_xmax - $PL_xmin) / $PL_nx;
  my $sizey = ($PL_ymax - $PL_ymin) / $PL_ny;
  my $x1 =  ($PL_xmin + ($x - 1) * $sizex);
  my $x2 = $x1 + $sizex - $PL_deltax;
  my $y1 = $PL_ymax - $y * $sizey;
  my $y2 = $y1 + $sizey - $PL_deltay;
  
  return($x1, $x2, $y1, $y2);
}

sub PanelCoordN
{
  my ($N) = @_;
  
  $N--;
  my $x = $N % $PL_nx + 1;
  my $y = int($N / $PL_nx) + 1;
  return PanelCoord($x, $y);
}
  

#########################################################

=item PlotPanelXY

  PlotPanelXY($win, $x, $y, %opt);

Plots a panel in the grid co-ordinated $x, $y. %opt can contain all possible options for xyplot.

=cut

sub PlotPanelXY
{
  my ($win, $x, $y, %opt) = @_;
  
  my ($x1, $x2, $y1, $y2);
  if(defined $opt{coord})
    {($x1, $x2, $y1, $y2) = @{$opt{coord}}}
  else
    {($x1, $x2, $y1, $y2) = PanelCoord($x, $y)}
  #print "$x1 $x2 $y1 $y2\n";die;
  
  
  $opt{BOX} ||= [0, 1, 0, 1];
  $opt{XBOX} = 'BST' unless defined $opt{XBOX};
  $opt{YBOX} = 'BST' unless defined $opt{YBOX};
  $opt{XLAB} ||= '';
  $opt{YLAB} ||= '';
  $opt{CHARSIZE} ||= 0.65;
  $opt{COLOR} ||= 'BLACK';
  $opt{MINTICKSIZE} = 0.5 unless defined $opt{MINTICKSIZE}; 
  $opt{MAJTICKSIZE} = 0.5 unless defined $opt{MAJTICKSIZE};

  $opt{VIEWPORT} = [$x1, $x2, $y1, $y2];
  $opt{XBOX} .= 'N' if $y == $PL_ny && !$opt{nonx} && !$opt{non};
  $opt{YBOX} .= 'N' if $x == 1 && !$opt{nony} && !$opt{non};
  $opt{XLAB} = '' if $y < $PL_ny && !$opt{forcexlab};
  $opt{YLAB} = '' if $x > 1 && !$opt{forceylab};
  if($opt{xbox})
  {
    $opt{XBOX} .= $opt{xbox};
    delete $opt{xbox}
  }
  if($opt{ybox})
  {
    $opt{YBOX} .= $opt{ybox};
    delete $opt{ybox}
  }
  
  
  delete $opt{forcexlab};
  delete $opt{forceylab};
  delete $opt{nonx};
  delete $opt{nony};
  delete $opt{non};
  delete $opt{coord};
  
  my $_x = pdl($opt{BOX}[0] - 1);
  my $_y = pdl($opt{BOX}[2] - 1);
  $win->xyplot($_x, $_y, %opt);
  
  return $x1, $x2, $y1, $y2
}

#########################################################

=item PlotPanelN

  PlotPanelN($win, $N, %opt);

=cut

sub PlotPanelN
{
  my ($win, $N, %opt) = @_;
  
  $N--;
  my $x = $N % $PL_nx + 1;
  my $y = int($N / $PL_nx) + 1;
  PlotPanelXY($win, $x, $y, %opt)
}

#########################################################

sub DistPlot
{
  my ($win, $data, %opt) = @_;
  
  my %defaults = (
    colour => [215, 215, 255],
    linewidth => 1,
    nbin => undef,
    min => undef,
    max => undef,
    width => 1,
    intbin => undef,
    symmetric => undef,
    boxes => undef,
    withmedian => undef,
    withtrimmedmean => undef,
    alpha => 0.5
  );
  %opt = parse(\%defaults, \%opt);
  my ($colour, $linewidth, $n, $min, $max, $width) = map {$opt{$_}} qw(colour  linewidth nbin min max width);
  
  if($opt{intbin})
  {
    ($min, $max) = minmax(pdl($data));
    $n = $max - $min + 1;
    $min -= 0.5;
    $max += 0.5;
  }
  
  my $num = @$data;
  
  my $R = zeroes(2) + $colour->[0];
  my $G = zeroes(2) + $colour->[1];
  my $B = zeroes(2) + $colour->[2];
  plscmap1($R, $G, $B);

  my $boxwidth = 0.8;

  for my $k (0 .. @$data - 1)
  {
    next unless defined $data->[$k];
    my $d = $data->[$k]->where(isgood($data->[$k]));
    $d->badflag(0);

    my ($y1, $y2, $h) = BuildHistogram2($d, n=>$n, min=>$min, max=>$max, renorm=>1);
    my $delta = at($y2,0) - at($y1,0);
    my $gap = 0.5 * (1 - $width) * $delta; 
    $y1 += $gap;
    $y2 -= $gap;
    
    my $scale = $boxwidth / max($h);
    my $N = nelem($y1);
    
    for my $i (0 .. $N - 1)
    {
      my $x = at($h, $i) * $scale;
      my ($x1, $x2);
      if($opt{symmetric})
      {
        $x2 = $k + $x/2;
        $x1 = $k - $x/2;
      }
      else
      {
        $x2 = $k + $boxwidth / 2;
        $x1 = $x2 - $x;
      }
      my $bottom = at($y1, $i);
      my $top = at($y2, $i);
      my $px = pdl($x2, $x2, $x1, $x1, $x2);
      my $py = pdl($bottom, $top, $top, $bottom, $bottom);
      
      #print "px = $px\npy = $py\n\n";
      
      plcol1(0);
      plfill($px, $py);
      plotboxelement($win, $px, $py, 'BLACK', 1, 1) if $opt{boxes};
    }
    #plotboxelement($win, pdl($x2, $x2), pdl(at($y1,0), at($y2,$N-1)), 'GREY', 1, 1);
    my $x0 = ($opt{symmetric}) ? $k : $k + $boxwidth / 2;
    LinePlot($win, pdl($x0, $x0),  pdl(at($y1,0), at($y2,$N-1)), COLOR=>'GREY');
    if($opt{withmedian})
    {
      my $med = median($d);
      my ($hm, $err) = interpolate($med, ($y1+$y2)/2, $h);
      my $hx = at($hm, 0) * $scale;
      LinePlot($win, pdl($k-$hx/2, $k+$hx/2), pdl($med, $med), COLOR=>'BLACK', LINEWIDTH=>3);
    }
    if($opt{withtrimmedmean})
    {
      my $tmean = TrimmedMean($d, $opt{alpha});
      my ($hm, $err) = interpolate($tmean, ($y1+$y2)/2, $h);
      my $hx = at($hm, 0) * $scale;
      LinePlot($win, pdl($k-$hx/2, $k+$hx/2), pdl($tmean, $tmean), COLOR=>'BLACK', LINEWIDTH=>3);
    }
  }
}

#########################################################

=item BoxPlot

  BoxPlot($win, $data, %opt);
  
Creates a boxplot (whisker plot) in window $win. $data is a reference to a list of piddles; each piddle cointains several data points (numbers) from which one individual box willl be created. There are following options available:

=over 4

=item B<colour>

Reference to a 3-element array, e.g. [0, 255, 255]. If not specified, each box will be white.

=item B<pcolour>

Colour of points in the cloud.

=item B<linewidth>

Line width.

=item B<cloud>

If defined (C<cloud => 1>), a cloud of points will be plotted on top of each box.

=item B<thickmed>

If defined, median is thicker than other lines.

=item B<boxwidth>

Relative width of each box. Default value is 0.8.

=back

B<Example>

  my $win = NewWindow();
  
  $win->xyplot(pdl(-10), pdl(-10), BOX => [-1, 2, -3, 3]);
  my $d1 = grandom(10);
  my $d2 = grandom(20) + 0.5;
  my $data = [$d1, $d2];
  BoxPlot($win, $data, colour => [0, 0, 255], linewidth => 2, could => 1);
  
  $win->close();

=cut

sub BoxPlot
{
  my ($win, $data, %opt) = @_;
  
  my %defaults = (
    colour => [255, 255, 255],
    pcolour => [0, 0, 0],
    lcolour => [0, 0, 0],
    linewidth => 1,
    thickmed => undef,
    psize => 0.7,
    pos => 0,
    solidwhisk => undef,
    cwidth => 0.1,
    cloud => undef,
    outliers => undef
  );
  %opt = parse(\%defaults, \%opt);
  my ($colour, $pcolour, $lcolour, $linewidth, $psize, $pos, $solidwhisk, $cwidth) = map {$opt{$_}} qw(colour  pcolour lcolour linewidth psize pos solidwhisk cwidth);
  
  my $medwidth = ($opt{thickmed}) ? 1 : 0;
  
  my $num = @$data;
  
  my $R = zeroes(2) + $colour->[0];
  my $G = zeroes(2) + $colour->[1];
  my $B = zeroes(2) + $colour->[2];
  plscmap1($R, $G, $B);

  my $boxwidth = ($opt{boxwidth}) ? $opt{boxwidth} : 0.8;
  for my $i (0 .. @$data - 1)
  {
    next unless defined $data->[$i];
    my $d = $data->[$i]->where(isgood($data->[$i]));
    $d->badflag(0);

    next if nelem($d) == 0;
    my ($min, $p05, $q1, $median, $q2, $p95, $max, $mean) = boxstats($d);
    
    my $ii = $pos + $i;
    my $x1 = $ii - $boxwidth / 2;
    my $x2 = $ii + $boxwidth / 2;
    my $sx1 = $ii - $boxwidth / 4;
    my $sx2 = $ii + $boxwidth / 4;

    my $x_ver = pdl($ii, $ii);
    my $x_hor = pdl($x1, $x2);
    my $x_shor = pdl($sx1, $sx2);
    my $x_box = pdl($x1, $x2, $x2, $x1, $x1);
    my $y_med = pdl($median, $median);
    my $y_box = pdl($q1, $q1, $q2, $q2, $q1);
    my $y_wh1 = pdl($p05, $q1);
    my $y_wh2 = pdl($q2, $p95);
    my $y_bot = pdl($p05, $p05);
    my $y_top = pdl($p95, $p95);

    #plscmap0($R, $G, $B);
    PDL::plcol1(0);
    unless($opt{nofill})
    {
      if($opt{swap})
        {PDL::plfill($y_box, $x_box)}
      else
        {PDL::plfill($x_box, $y_box)}
    }

    if($opt{cloud})
    {
      #my $xcloud = 0.5 * (random($data->[$i]) - 0.5) + $i;
      my $xcloud = $cwidth * grandom($d) + $ii;
      #print $data->[$i]->(0:10), $xcloud(0:10), "\n";
      my ($x_, $y_) = ($opt{swap}) ? ($d, $xcloud) : ($xcloud, $d);
      $win->xyplot($x_, $y_,
        PLOTTYPE => 'POINTS',
        XBOX => '', YBOX => '', XLAB => '', YLAB => '',
        SYMBOL => 751,
        SYMBOLSIZE => $psize,
        COLOR => $pcolour,
      );  
    }
    
    if($opt{outliers})
    {
      my $d = $data->[$i];
      my $yout = $d->where(($d <= $p05) | ($d >= $p95));
      my $xout =  zeroes($yout) + $ii;
      my ($x_, $y_) = ($opt{swap}) ? ($yout, $xout) : ($xout, $yout);
      $win->xyplot($x_, $y_,
        PLOTTYPE => 'POINTS',
        XBOX => '', YBOX => '', XLAB => '', YLAB => '',
        SYMBOL => 751,
        SYMBOLSIZE => $psize,
        COLOR => $pcolour,
      );  
    }


    #my $colour = 'BLACK';
    my $s = $opt{swap};
    my $sw = ($solidwhisk) ? 1 : 4;
    plotboxelement($win, $x_hor, $y_med, $lcolour, $linewidth+$medwidth, 1, $s);
    plotboxelement($win, $x_box, $y_box, $lcolour, $linewidth, 1, $s);
    plotboxelement($win, $x_ver, $y_wh1, $lcolour, $linewidth, $sw, $s);
    plotboxelement($win, $x_ver, $y_wh2, $lcolour, $linewidth, $sw, $s);
    plotboxelement($win, $x_shor, $y_bot, $lcolour, $linewidth, 1, $s);
    plotboxelement($win, $x_shor, $y_top, $lcolour, $linewidth, 1, $s);
   
   
    if($opt{mean})
    {
      my ($x_, $y_) = ($opt{swap}) ? (pdl($mean), pdl($ii)) : (pdl($ii), pdl($mean));
      $win->xyplot($x_, $y_,
        PLOTTYPE => 'POINTS',
        XBOX => '', YBOX => '', XLAB => '', YLAB => '',
        SYMBOL => 751,
        SYMBOLSIZE => 1.1,
        COLOR => $pcolour,
      );  
    }
  }
}

sub plotboxelement
{
  my ($win, $x, $y, $colour, $linewidth, $linestyle, $swap) = @_;
  
  my ($x_, $y_) = ($swap) ? ($y, $x) : ($x, $y);
  
  $win->xyplot($x_, $y_,
      PLOTTYPE => 'LINE',
      XBOX => '', YBOX => '', XLAB => '', YLAB => '',
      COLOR => $colour,
      LINEWIDTH => $linewidth,
      LINESTYLE => $linestyle
  );  
}

sub boxstats
{
  my $x = shift;
  
  my $n = nelem($x);
  my $xs = qsort $x;
  my $min = at($xs, 0);
  my $p05 = percentile($xs, 0.05);
  my $q1 = percentile($xs, 0.25);
  my $median = percentile($xs, 0.5);
  my ($mean) = stats($x);
  my $q2 = percentile($xs, 0.75);
  my $p95 = percentile($xs, 0.95);
  my $max = at($xs, $n - 1);
  
  return ($min, $p05, $q1, $median, $q2, $p95, $max, $mean);
}

sub percentile
{
  my ($x, $p) = @_;
  
  my $n = nelem($x);
  my $i = $p * ($n - 1);
  my $int = int($i);
  my $per;
  if($i == $int)
    {$per = at($x, $i)}
  elsif($i < $n - 1)
    {$per = (at($x, $i) + at($x, $i+1)) / 2}
  else
    {$per = at($x, $n-1)}
  return $per;
}

#######################################################

=item BuildHistogram2

  my ($x1, $x2, $h) = BuildHistogram($d, %opt);

Options: n, min, max, renorm, extend.

Builds a histogram from data $d with $n bars, between $min and $max (optional), renormalized to its integral (if $renorm = 1) or sum ($renorm = 2). Returns left and right bar limits, $x1 and $x2, and the histogram values $h.

=cut

sub BuildHistogram2
{
  my ($d, %opt) = @_;
  
  my %defaults = (
    n => undef,
    min => undef,
    max => undef,
    renorm => undef,
    extend=>undef,
  ); 
  %opt = parse(\%defaults, \%opt);
  my ($n, $min, $max, $renorm, $extend) = map {$opt{$_}} qw(n min max renorm extend);
  
  my ($mean, $stdev, $med, $dmin, $dmax) = stats($d);
  unless(defined $n)
  {
    my $num = nelem($d);
    my $h = 3.49 * $stdev / $num**(1/3); # Scott's formula for suggested bin witdth
    $n = int(($dmax - $dmin) / $h) + 1;
  } 

  ($min, $max) = ($dmin, $dmax) unless defined $min;
  my $step = ($max - $min) / $n;
  my $x1 = $min + sequence($n) * $step;
  my $x2 = $x1 + $step;

  my $h = histogram($d, $step, $min, $n);
  if($renorm)
  {
    $h /= $step * sum($h) if $renorm == 1;
    $h /= sum($h) if $renorm == 2;
  }
  
  return ($x1, $x2, $h);
}

#######################################################

=item BuildHistogram

  my ($x, $h) = BuildHistogram($d, %opt);
  
Options: n, min, max, renorm, extend.
  
Builds a histogram from data $d with n bars, between min and max, renormalized to its integral (if renorm = 1) or sum (renorm = 2). Returns middle of bar x-positons, $x, and the histogram values $h.

If n is not given, it will be calculated using Scott's formula.

=cut

sub BuildHistogram
{
  my ($d, %opt) = @_;
  
  
  my %defaults = (
    n => undef,
    min => undef,
    max => undef,
    renorm => undef,
    extend=>undef,
  ); 
  %opt = parse(\%defaults, \%opt);
  my ($n, $min, $max, $renorm, $extend) = map {$opt{$_}} qw(n min max renorm extend);
  
  my ($mean, $stdev, $med, $dmin, $dmax) = stats($d);
  
  unless(defined $n)
  {
    my $num = nelem($d);
    my $h = 3.49 * $stdev / $num**(1/3); # Scott's formula for suggested bin witdth
    $n = int(($dmax - $dmin) / $h) + 1;
    if($extend)
    {
      $n += 2;
      $dmin -= $h;
      $dmax += $h;
    }
  } 

  ($min, $max) = ($dmin, $dmax) unless defined $min;
  my $step = ($max - $min) / $n;
  my $x = $min + sequence($n) * $step + 0.5 * $step;

  my $h = histogram($d, $step, $min, $n);
  #print "HISTOGRAM: step = $step, min = $min, n = $n\n";
  my $norm = 1;
  if($renorm)
  {
    $norm = $step * sum($h) if $renorm == 1;
    $norm = sum($h) if $renorm == 2;
    $h /= $norm;
  }

  return ($x, $h, $norm);
}

#########################################################

=item KDE

  my ($x, $y) = KDE($d, $h, $size, $N);
  
Builds kernel density estimator from sample $d, with window width $h, $size standard deviations around the mean, using $N points. Gaussian kernel is used.
 
=cut

sub KDE
{
  my ($d, $h, $size, $N) = @_;
  
  my ($m, $s) = stats($d);
  my $n = nelem($d);
  $h ||= 1.06 * $s * $n**-0.2;  # good for normal distribution
  $size ||= 5;
  $N ||= 300;
  #print "Using h = $h\n";
  
  my $x = $size * $s * (2*sequence($N)/($N - 1) - 1) + $m;
  my @y = ();
  for my $k (0 .. $N - 1)
  {
    my $u = (at($x, $k) - $d) / $h;
    my $f = 1 / sqrt(2*$PI) * sum(exp(-0.5*$u*$u));
    push @y, $f / ($n * $h);
  }
  my $y = pdl(@y);
  return $x, $y; 
}

#########################################################


sub MultiBarPlot
{
  my ($win, $x, $y, $cols, %opt) = @_;
  
  my $width = $opt{width};
  $width ||= 1;
  
  my $gap = (1 - $width) * (at($x,1) - at($x,0));  # equidistant points
  my $n = nelem($x);
  
  my $l = 0.5 * ($x(0:$n-2) + $x(1:$n-1) + $gap);
  my $r = 0.5 * ($x(0:$n-2) + $x(1:$n-1) - $gap);
  
  my $first = 2 * at($x,0) - at($r,0);
  my $last =  2 * at($x,$n-1) - at($l,$n-2);
  
  my $x1 = append(pdl($first), $l);
  my $x2 = append($r, pdl($last));
  
  my @xmid = ();
  my $N = @$y;
  my $delta = (at($x2, 0) - at($x1, 0)) / $N;
  for my $i (0 .. $N - 1)
  {
    my $xx1 = $x1 + $i * $delta;
    my $xx2 = $x1 + ($i + 1) * $delta;
    push @xmid, ($xx1 + $xx2) / 2;
    #print "$i, $delta, $xx1, $xx2, $y->[$i]\n";
    BarPlot2($win, $xx1, $xx2, $y->[$i], %opt, colour => $cols->[$i]);
  }
  return \@xmid;
}

#########################################################

=item BarPlot

  BarPlot($win, $x, $y, %opt);
  
Creates a bar plot in window $win. $x and $y are histogram data (equidistant x points assumed). Options are the same as in BarPlot2, plus

=over 4

=item B<width>

Fractional width of each bar. If not defined, a default value of 1 is assumed, i.e. there is no gap between the bars.

=back

=cut

sub BarPlot
{
  my ($win, $x, $y, %opt) = @_;
  
  my $width = $opt{width};
  $width ||= 1;
  
  my $gap = (1 - $width) * (at($x,1) - at($x,0));  # equidistant points
  my $n = nelem($x);
  
  my $l = 0.5 * ($x(0:$n-2) + $x(1:$n-1) + $gap);
  my $r = 0.5 * ($x(0:$n-2) + $x(1:$n-1) - $gap);
  
  my $first = 2 * at($x,0) - at($r,0);
  my $last =  2 * at($x,$n-1) - at($l,$n-2);
  
  my $x1 = append(pdl($first), $l);
  my $x2 = append($r, pdl($last));
  
  BarPlot2($win, $x1, $x2, $y, %opt);
}

#########################################################

=item BarPlot2

  BarPlot2($win, $x1, $x2, $y, %opt);
  
Creates a bar plot of a histogram created by C<BuildHistogram> (or by other means), in window $win. $x1 and $x2 are left and right ends of bars, $y is hights of bars. There are the following options available:

=over 4

=item B<colour>

Reference to a 3-element list [R, G, B]. If not defined, [100, 100, 100] is used.

=item B<boxes>

If defined (C<boxes => 1>), a black box is drawn around each bar.

=item B<y0>

The bottom y-coordinate of the bars.

=back

=cut

sub BarPlot2
{
  my ($win, $x1, $x2, $y, %opt) = @_;
  
  my $colour = $opt{colour};
  my $colours = $opt{colours};
  my $boxlinewidth = ($opt{boxlinewidth}) ? $opt{boxlinewidth} : 1;
  my $boxlinecolour = ($opt{boxlinecolour}) ? $opt{boxlinecolour} : 'BLACK';
  if(defined $colour && ref($colour) !~ /ARRAY/)
  {
    die "Colour $colour unrecognized" unless exists $colours{$colour};
    $colour = $colours{$colour}
  }
  $colour ||= [100, 100, 100];

  my $y0 = $opt{y0};
  $y0 = -10 unless defined $y0;

  my $R = pdl($colour->[0], 0);
  my $G = pdl($colour->[1], 0);
  my $B = pdl($colour->[2], 0);
  plscmap1($R, $G, $B);

  my $N = nelem($x1);
  for my $i (0 .. $N - 1)
  {
    if(defined $colours)
    {
      my @col = @{$colours->[$i]};
      my $R = pdl($col[0], 0);
      my $G = pdl($col[1], 0);
      my $B = pdl($col[2], 0);
      plscmap1($R, $G, $B);
    }
    
    my $left = at($x1, $i);
    my $right = at($x2, $i);
    my $hx = pdl($left, $right, $right, $left, $left);
    #my $py = pdl($y0, $y0, $y0 + at($y,$i), $y0 + at($y,$i), $y0);
    my $hy = pdl($y0, $y0, at($y,$i), at($y,$i), $y0);
    my ($px, $py) = ($opt{swap}) ? ($hy, $hx) : ($hx, $hy);
    plcol1(0);
    plfill($px, $py);
    plotboxelement($win, $px, $py, $boxlinecolour, $boxlinewidth, 1) if $opt{boxes};
  }
}

#########################################################

sub HistPlot
{
  my ($win, $x, $y, %opt) = @_;
  
  my $n = nelem($x);
  
  my $l = 0.5 * ($x(0:$n-2) + $x(1:$n-1));
  my $r = 0.5 * ($x(0:$n-2) + $x(1:$n-1));
  
  my $first = 2 * at($x,0) - at($r,0);
  my $last =  2 * at($x,$n-1) - at($l,$n-2);
  
  my $x1 = append(pdl($first), $l);
  my $x2 = append($r, pdl($last));
  
  HistPlot2($win, $x1, $x2, $y, %opt);
}

#########################################################

=item HistPlot

  HistPlot($win, $x1, $x2, $y, %opt);
  
Creates a line plot of a histogram created by C<BuildHistogram> (or by other means) in window $win. $x1 and $x2 are left and right ends of bars, $y is hights of bars. There are the following options available:

=over 4

=item B<colour>

Reference to a 3-element list [R, G, B] or a text name (e.g. 'BLACK').

=item B<linewidth>

Line width.

=item B<linestyle>

Line style.

=back

=cut

sub HistPlot2
{
  my ($win, $x1, $x2, $y, %opt) = @_;
  
  my $colour = $opt{colour};
  my $linewidth = $opt{linewidth};
  my $linestyle = $opt{linestyle};
  
  $colour ||= 'BLACK';
  $linewidth ||= 2;
  $linestyle ||= 1;

  my $N = nelem($y);
  my $M = 3 * $N;
  my $hx = zeroes($M);
  my $hy = zeroes($M);

  $hx(0:$M-1:3) .= $x1;
  $hx(1:$M-1:3) .= $x1;
  $hx(2:$M-1:3) .= $x2;

  $hy(0) .= $y(0);
  $hy(3:$M-1:3) .= $y(0:$N-2);
  $hy(1:$M-1:3) .= $y;
  $hy(2:$M-1:3) .= $y;

  my ($px, $py) = ($opt{swap}) ? ($hy, $hx) : ($hx, $hy);

  $win->xyplot($px, $py,
      PLOTTYPE => 'LINE', LINESTYLE => 1,
      COLOR => $colour, LINEWIDTH => $linewidth,
      XBOX => '', YBOX => '', XLAB => '', YLAB => '',
  );
}

############################

=item BoxMinMax

  my ($boxmin, $boxmax) = BoxMinMax($x, $spacer);
  
=cut

sub BoxMinMax
{
  my ($x, $spacer) = @_;
  $spacer = 0.1 if !defined $spacer;
  
  my ($min, $max) = minmax($x);
  my $delta = $spacer * ($max - $min);
  my $boxmin = $min - $delta;
  my $boxmax = $max + $delta;
  return($boxmin, $boxmax);
}

############################

=item GraphColour

  my $colour = GraphColour($i);
  
=cut

sub GraphColour
{
  my $i = shift;
  my @col = qw(BLACK RED BLUE GREEN TURQUOISE AQUAMARINE PINK MAGENTA WHEAT GREY BROWN);
  
  return $col[$i % $#col]
}

############################

=item OverPlot

  OverPlot($win, $x, $y, %opt);
  
=cut

sub OverPlot
{
  my ($win, $x, $y, %opt) = @_;
  
  $win->xyplot($x, $y, %opt,
    XBOX => '',
    YBOX => '',
    XLAB => '',
    YLAB => ''
  );
}

############################

=item PointPlot

  PointPlot($win, $x, $y, %opt);
  
=cut

sub PointPlot
{
  my ($win, $x, $y, %opt) = @_;
  
  $opt{SYMBOLSIZE} ||= 0.8;
  $opt{COLOR} ||= 'BLACK';
  $opt{SYMBOL} ||= 751;

  $win->xyplot($x, $y, %opt,
    PLOTTYPE => 'POINTS',
    XBOX => '',
    YBOX => '',
    XLAB => '',
    YLAB => ''
  );
}

############################

=item LinePlot

  LinePlot($win, $x, $y, %opt);
  
=cut

sub LinePlot
{
  my ($win, $x, $y, %opt) = @_;
  
  $opt{LINEWIDTH} ||= 1;
  $opt{LINESTYLE} ||= 1;
  $opt{COLOR} ||= 'BLACK';

  $win->xyplot($x, $y, %opt,
    PLOTTYPE => 'LINE',
    XBOX => '',
    YBOX => '',
    XLAB => '',
    YLAB => ''
  );
}

############################

=item 

  ShadeErrPlot($win, $x, $y, $sy, %opt);
  
Plot shaded error bar area around. Should be used before a line is plotted. Options are:

  colour - shade colour
  fade - fading factor: 0 - no fade, 1 - full fade. Default value is 0.2
  
Usage:

  ShadeErrPlot($w, $x, $y, $sy, color=>$col, fade=>0.3); 

=cut

sub ShadeErrPlot
{
  my ($w, $x, $y, $sy, %opt) = @_;
  
  my $col = $opt{colour};
  my $fade = $opt{fade};
  
  $col ||= 'BLACK';
  $fade ||= 0.2;
  
  my $px = append($x, $x(-1:0));
  my $py = append($y - $sy, $y(-1:0) + $sy(-1:0));
  
  my $shade;
  $shade->[$_] = 255 - (255 - $col->[$_]) * $fade for (0 .. 2);
  plscmap1(pdl($shade->[0],0), pdl($shade->[1],0), pdl($shade->[2],0));
  plcol1(0);
  plfill($px, $py);
}


############################

=item XYShadePlot

  XYShadePlot($w, $x, $y, %opt);
  
  Plot a shaded scatter plot. Options are
  - xmin, xmax, ymin, ymax - define the box
  - xnbin, ynbin, nbin - nuber of bins in x, y or both
  - xextend, yextend - extension area for histogram to avoid piling up at edges
  - gamma - power index for scaling law (good numbers are 0.1-0.3)
  - hue - hue of the plot (0 - 255)

=cut

sub XYShadePlot
{
  my ($w, $x, $y, %opt) = @_;
  
  my %defaults = (
    xmin => undef,
    xmax => undef,
    ymin => undef,
    ymax => undef,
    xextend=>0,
    yextend => 0,
    xnbin => 100,
    ynbin => 100,
    nbin => undef,
    gamma => 0.1,
    hue => 230,
    minlum => 0
  ); 
  %opt = parse(\%defaults, \%opt);
  my ($xmin, $xmax, $ymin, $ymax, $nbin, $xnbin, $ynbin, $xextend, $yextend, $gamma, $hue, $minlum) = map {$opt{$_}} qw(xmin xmax ymin ymax nbin xnbin ynbin xextend yextend gamma hue minlum);

  if(defined $nbin)
  {
    $xnbin = $nbin;
    $ynbin = $nbin
  }

  my $xstep = ($xmax - $xmin + 2*$xextend) / $xnbin;
  my $ystep = ($ymax - $ymin + 2*$yextend) / $ynbin;
  
  my $h = histogram2d($x, $y, $xstep, $xmin-$xextend, $xnbin, $ystep, $ymin-$yextend, $ynbin);
  $h /= max($h);
  $h = $h**$gamma;  # rescale for good colour scale  
  
  my $pos = pdl(0, 1);
  my $H = pdl($hue, $hue);
  my $L = pdl(1, $minlum);
  my $S = pdl(1,  1);
  PDL::plscmap1l(0, $pos, $H, $L, $S, pdl[]);
  #$w->shadeplot ($h, 128, BOX => [$xmin, $xmax, $ymin, $ymax], ZRANGE => [0,1]);
  
  HeatMapRows($h, {xmin=>$xmin-$xextend, xmax=>$xmax+$xextend, ymin=>$ymin-$yextend, ymax=>$ymax+$yextend});
  
}

############################

=item MAPlot

  MAPlot($w, $pan, $x, $y, %opt);
  
  Plot a shaded MA plot.
  - xlab, ylab - name for x and y data
  - xmin, xmax, ymin, ymax - define the box
  - xnbin, ynbin, nbin - nuber of bins in x, y or both
  - xextend, yextend - extension area for histogram to avoid piling up at edges
  - gamma - power index for scaling law (good numbers are 0.1-0.3)
  - hue - hue of the plot (0 - 255)

=cut

sub MAPlot
{
  my ($w, $pan, $x, $y, %opt) = @_;
  
  my %defaults = (
    plottype => 'shade',
    xmin => undef,
    xmax => undef,
    ymin => undef,
    ymax => undef,
    xlab => 'x',
    ylab => 'y',
    perc => 0.95,
    symbolsize => 0.1,
    charsize => 0.6,
    nbin=>160,
    xextend=>undef,
    yextend=>undef,
    gamma=>0.3,
    hue=>250
  ); 
  %opt = parse(\%defaults, \%opt);
  my ($plottype, $xmin, $xmax, $ymin, $ymax, $xlab, $ylab, $symbolsize, $nbin, $xextend, $yextend, $gamma, $hue, $perc, $charsize) = map {$opt{$_}} qw(plottype xmin xmax ymin ymax xlab ylab symbolsize nbin xextend yextend gamma hue perc charsize);
  
  warn "Negative data in MAPlot; ignored\n" if any(($x <0) & ($y < 0));
  my $good = which(isgood($x) & isgood($y) & ($x > 0) & ($y > 0));
  $x = $x($good);
  $y = $y($good);
  
  my $corr = corr($x, $y);
  my $mean = log10(($x + $y) / 2);
  my $rat = log($y / $x) / log(2);
  my $med = median($rat);
  my $alpha = (1 - $perc) / 2;
  my $lo = pct($rat, $alpha);
  my $up = pct($rat, 1 - $alpha);
  
  #printf "lo = %.2g   med = %.2g   up = %.2g\n", $lo, $med, $up;
  #my $F = (2**($up) + 2**(-$lo))/2;
  my $F = ($up - $lo)/2;
  $F = sqrt(avg($rat**2));
  
  my ($_xmin, $_xmax) = BoxMinMax($mean);
  my ($_ymin, $_ymax) = BoxMinMax($rat);

  $xmin = $_xmin unless defined $xmin;
  $xmax = $_xmax unless defined $xmax;  
  $ymin = $_ymin unless defined $ymin;
  $ymax = $_ymax unless defined $ymax;
  
  $xextend = ($xmax - $xmin) / $nbin unless defined $xextend;
  $yextend = ($ymax - $ymin) / $nbin unless defined $yextend;
  
  PlotPanelN($w, $pan, BOX => [$xmin, $xmax, $ymin, $ymax], %PL_empty);
  
  if($plottype eq 'point')
    {PointPlot($w, $mean, $rat, SYMBOLSIZE=>$symbolsize)}
  elsif($plottype eq 'shade')
    {XYShadePlot($w, $mean, $rat, xmin=>$xmin, xmax=>$xmax, ymin=>$ymin, ymax=>$ymax, nbin=>$nbin, xextend=>$xextend, yextend=>$yextend, gamma=>$gamma, hue=>$hue)}
  else
    {carp "Huh? Unknown plottype=$plottype in MAPlot\n"}
  #LinePlot($w, $xx, $yy, COLOR=>'RED');

  LinePlot($w, pdl($xmin, $xmax), pdl($med, $med), COLOR=>'RED', LINEWIDTH=>3);
  LinePlot($w, pdl($xmin, $xmax), pdl($lo, $lo), COLOR=>'RED', LINEWIDTH=>1);
  LinePlot($w, pdl($xmin, $xmax), pdl($up, $up), COLOR=>'RED', LINEWIDTH=>1);

  PlotPanelN($w, $pan,
    BOX => [$xmin, $xmax, $ymin, $ymax],
    XLAB => "log#d10#u ($xlab+$ylab)/2",
    YLAB => "log#d2#u $ylab/$xlab",
    forcexlab=>1,
    forceylab=>1,
    non=>1,
    XBOX => 'BINSTF', YBOX => 'BINSTF',
    CHARSIZE => $charsize
  );
  
  return $F, $corr;
}


############################

=item TopText

  TopText($win, $s, %opt);
  
=cut

sub TopText
{
  my ($win, $s, %opt) = @_;
  
  $opt{CHARSIZE} ||= 0.8;
  $opt{COLOR} ||= 'BLACK';
  
  my $x = (defined $opt{x}) ? $opt{x} : 0.05;
  my $y = (defined $opt{y}) ? $opt{y} : -1;
  my $c = (defined $opt{c}) ? $opt{c} : 0;
  my $pos = (defined $opt{bottom}) ? 'b' : 't';
  
  delete $opt{x};
  delete $opt{y};
  delete $opt{c};
  delete $opt{bottom};

  $win->text($s, %opt, TEXTPOSITION=>[$pos, $y, $x, $c]);
}


=item FullBox

  FullBox($win);
  
=cut

sub FullBox
{
  my $w = shift;
  $w->xyplot(pdl(-1),pdl(-1), VIEWPORT=>[0,1,0,1], BOX=>[0,1,0,1], XBOX=>'', YBOX=>'', XLAB=>'', YLAB=>'');
}

##########################################################

=item PLStore

  my $PL = PLStore();
  
  Stores all global PL parameters ($PL_nx, $PL_xmin and so on) in one variable, so they can be played with and then restored.
  
=cut

sub PLStore
{
  my $PL = [$PL_nx, $PL_ny, $PL_xmin, $PL_xmax, $PL_ymin, $PL_ymax, $PL_deltax, $PL_deltay];
  
  return $PL
}

=item PLRecover

  PLRecover($PL);
  
  Recovers all global PL parameters stored with PLStore().
  
=cut

sub PLRecover
{
  my $PL = shift;
  
  ($PL_nx, $PL_ny, $PL_xmin, $PL_xmax, $PL_ymin, $PL_ymax, $PL_deltax, $PL_deltay) = @$PL; 
}

1;
