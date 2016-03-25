package Heatmap;
require Exporter;
require PDL;
require PDL::NiceSlice;

=head1 NAME

B<Heatmap> - heat maps for PDL and PLplot

=head1 SYNPOSIS

  use PDL;
  use PDL::Graphics::PLplot;
  use CompBio::Heatmap;

  my $win = PDL::Graphics::PLplot->new(DEV => 'xwin');
  $win->xyplot(pdl(0,10), pdl(0,100));
  
  SCMapGreenBlackRed();
  HeatMapRows($map);
  HeatMapRows($map, {xmax => $nx, ymax => $ny, renorm => 1});
  
=head1 DESCRIPTION

Plots a heatmap with an arbitrary colour map in a viewport created by, e.g. xyplot. The input map is a two-dimensional piddle, with values normalized between 0 and 1. Values outside [0, 1] will cause errors.
  
=head1 FUNCTIONS

=over 4
  
=cut

use strict;
use warnings;

use PDL;
use PDL::NiceSlice;
use Switch 'Perl6';

our $VERSION = 0.31;

our @ISA = ("Exporter");
our @EXPORT = qw(
  HeatMapRows
  SCMap
  SCMapGreenBlackRed
  SCMapGreyscale
  SCMapReverseGreyscale
  SCMapHeat
  SCHuescaleMap
  plscmap1r
);


########################################################

sub HeatMapRows
{
  my ($map, $opt) = @_;

  my ($ncols, $nrows) = dims $map;
  
  my $xmin = $opt->{xmin} || 0;
  my $xmax = (defined $opt->{xmax}) ? $opt->{xmax} : $ncols - 1;
  my $ymin = $opt->{ymin} || 0;
  my $ymax = (defined $opt->{ymax}) ? $opt->{ymax} : $nrows - 1;
  my $ydel = ($ymax - $ymin) / $nrows;
  
  # $opt->{renorm} means renormalize each profile independelntly
  # here we renormalize all of them together
  if(!defined $opt->{renorm} && !defined $opt->{min} && !defined $opt->{max})
  {
    my ($min, $max) = minmax($map);
    $map = ($map - $min) / ($max - $min) if $max > $min;
  }
  
  if(defined $opt->{min} && defined $opt->{max})
  {
    $map = ($map - $opt->{min}) / ($opt->{max} - $opt->{min}) if $opt->{max} > $opt->{min};
  }
  
  for my $i (0 .. $nrows - 1)
  {
    my $p = $map(,($i));   # $i-th row ($i-th profile)
    my $y = $ymin + $i * $ydel;
    OneHeatRow($p, $y, $ydel, $xmin, $xmax, $opt);
  }
}

########################################################

sub OneHeatRow
{
  my ($p, $y, $dy, $xmin, $xmax, $opt) = @_;
  
  if($opt->{renorm})
  {
    my ($min, $max) = minmax($p);
    $p = ($p - $min) / ($max - $min) if $max > $min;
  }
  else
  {
    $p->where($p > 1) .= 1;
    $p->where($p < 0) .= 0;
  }
  my $N = nelem($p);
  $xmin = 0 unless defined $xmin;
  $xmax = $N - 1 unless defined $xmax;
  #print "===$p\n";
  my $delta = ($xmax - $xmin) / $N;
  #print "xmin=$xmin xmax=$xmax  delta=$delta  N=$N\n";
  my $py = pdl($y, $y, $y + $dy, $y + $dy);
  for my $i (0 .. $N - 1)
  {
    my $x = $xmin + $delta * $i;
    my $px = pdl($x, $x + $delta, $x + $delta, $x);
    PDL::plcol1(at($p, $i));
    PDL::plfill($px, $py);
  }
}

########################################################

sub SCMapGreenBlackRed
{
  my $pos = pdl(0, 0.5, 1);
  my $R = pdl(0, 0, 1);
  my $G = pdl(1, 0, 0);
  my $B = pdl(0, 0, 0);
  PDL::plscmap1l(1, $pos, $R, $G, $B, pdl[]);
}  

########################################################

sub SCMapGreyscale
{
  my $pos = pdl(0, 1);
  my $R = pdl(0, 1);
  my $G = pdl(0, 1);
  my $B = pdl(0, 1);
  PDL::plscmap1l(1, $pos, $R, $G, $B, pdl[]);
}  
  
########################################################

sub SCMapReverseGreyscale
{
  my $pos = pdl(0, 1);
  my $R = pdl(1, 0);
  my $G = pdl(1, 0);
  my $B = pdl(1, 0);
  PDL::plscmap1l(1, $pos, $R, $G, $B, pdl[]);
}  
  
########################################################

sub SCMapHeat
{
  my $pos = pdl(0, 0.5, 1);
  my $R = pdl(0, 1, 1);
  my $G = pdl(0, 0, 1);
  my $B = pdl(0, 0, 0);
  PDL::plscmap1l(1, $pos, $R, $G, $B, pdl[]);
}

########################################################

sub SCHuescaleMap
{
  my ($hue, $minlum, $reverse) = @_;
  
  $hue = 230 unless defined $hue;
  $minlum = 0.15 unless defined $minlum;
  my $pos = pdl(0, 1);
  my $H = pdl($hue, $hue);
  my $L = ($reverse) ? pdl(1, $minlum) : pdl($minlum, 1);
  my $S = pdl(1, 1);
  PDL::plscmap1l(0, $pos, $H, $L, $S, pdl[]);
}  

########################################################

sub SCMap
{
  my $hmap = shift;
  given($hmap)
  {
    when 'gbr'    {SCMapGreenBlackRed()}
    when 'heat'   {SCMapHeat()}
    when 'grey'   {SCMapGreyscale()}
    when 'rgrey'  {SCMapReverseGreyscale()}
    when 'blue'   {SCMapBluescale()}
    when 'rblue'  {SCMapReverseBluescale()}
  }
}

########################################################

=item plscmap1r

  plscmap1r($pos, $r, $g, $b);
  
An alternative to standard function plscmap1l where interpolation is in RGB space, as opposed to hue space.

=cut

sub plscmap1r
{
  my ($pos, $r, $g, $b) = @_;
  
  my $n = nelem($pos);
  die "Need n > 1\n" unless $n > 1;
  die "Wrong number of elements in red\n" unless nelem($r) == $n;
  die "Wrong number of elements in green\n" unless nelem($g) == $n;
  die "Wrong number of elements in blue\n" unless nelem($b) == $n;
  
  my $N = 255;
  
  my $x = sequence($N) / ($N - 1);
  my ($R) = interpolate($x, $pos, $r);
  my ($G) = interpolate($x, $pos, $g);
  my ($B) = interpolate($x, $pos, $b);
  
  PDL::plscmap1(255*$R, 255*$G, 255*$B);
}

1;


  

