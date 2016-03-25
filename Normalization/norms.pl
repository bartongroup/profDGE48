#!/sw/bin/perl

=head1 NAME

B<norms.pl>

=head1 DESCRIPTION

Compares various normalization methods. Also, creates files with normalizing factors that are used by C<plot_gene_pileup.pl> and other scripts. There are no options.

=cut


use strict;
use warnings;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use GRNASeq;
use GetOpt;

use Stats;
use PLgraphs;
use Tools;

$| = 1;

Stamp();
ReadDefs();
GetMyOptions();

#my @norms = qw(deseq tmm totcount totcountm);

my @norms = qw(deseq totcountm);
my @cols = qw(BLUE GREEN RED GOLD2 GREY);

my %h = ();
my %x = ();
for my $norm (@norms)
{
  #my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel, $f, $xsel) = 
  my %d = GetNormalizedData(
    norm=>$norm,
    nonzero=>1,
    type=>$type,
    clean=>$clean,
    spclean=>$spclean
  );
  for my $cond (@conds)
  {
    my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel, $f, $xsel) = @{$d{$cond}};
    $h{$cond}{$norm} = $f;
    $x{$cond} = $xsel + 1;
  }
}


PlotNorms();
CompareNorms();


sub PlotNorms
{
  my $w = NewWindow($psfile);
  $PL_ny = 2;
  
  my $pan = 1;
  my $nrep = 48;
  my $maxx = $nrep + 1;
  for my $cond (@conds)
  {
    my $x = $x{$cond};
    my ($min, $max) = (0, 1.99);
    PlotPanelN($w, $pan++,
      BOX => [0, $maxx, $min, $max],
      XLAB => 'Replicate #',
      YLAB => 'Normalization factor',
      #xbox => 'I', ybox => 'I'
    );
    
  #  for my $r (1 .. $nrep)
  #    {LinePlot($w, pdl($r,$r), pdl($min, $max), COLOR=>'GREY')} 
    
    my $i = 0;
    for my $norm (@norms)
    {
      my $f = $h{$cond}{$norm};  
      my $col = $cols[$i];
      
      print nelem($x), " ",nelem($f), "\n";
      
      LinePlot($w, $x, $f, COLOR=>'GREY');
      PointPlot($w, $x, $f, COLOR=>$col, SYMBOLSIZE=>0.9);
    
      wcols "%2d %5.3f", $x, $f, "${cond}_${norm}_normfac.txt";
    
      my $y = 0.5 - 0.08 * $i;
      LinePlot($w, pdl(2,5), pdl($y, $y), COLOR=>$col);
      $w->text($norm, COLOR=>$col, CHARSIZE=>0.6, TEXTPOSITION=>[5.5, $y, 0, 0, 0]);
      $i++;
    }
    
    LinePlot($w, pdl(0,$maxx+1), pdl(1,1), COLOR=>'RED');
    TopText($w, $cond, x=>0.07, y=>-2, CHARSIZE => 0.7);
    
    my $spf = SpikeinFile('yeast', $cond);
    #my ($rep, $sp, $spe) = rcols $spf, 0, 1, 2;
    #my $sp_norm = 1 / $sp;
    #my $sp_err = $spe / $sp**2;
    my $sp = rcols $spf, 1;
    my $sp_norm = 1 / $sp;
    
    #PointPlot($w, $x, $sp_norm($x-1), YERRORBAR=>2*$sp_err($x-1));
    PointPlot($w, $x, $sp_norm($x-1));
  }

  $w->close();
}

sub CompareNorms
{
  my $w = NewWindow($psfile);
  
  $PL_nx = 2;
  $PL_ymin = 0.5;
  
  my $pan = 1;
  my $norm1 = 'totcountm';
  my $norm2 = 'spikein';
    
  for my $cond (@conds)
  {
    my ($min, $max) = (0, 3.5);
    PlotPanelN($w, $pan++,
      BOX => [$min, $max, $min, $max],
      XLAB => $norm1,
      YLAB => $norm2,
      #xbox => 'I', ybox => 'I'
    );
    
    my $sel = $x{$cond} - 1;
    
    
    my $x = $h{$cond}{$norm1};
    my ($y, $ye);
    
    if($norm2 eq 'spikein')
    {
      my $spf = SpikeinFile('yeast', $cond);
      my ($rep, $sp, $spe) = rcols $spf, 0, 1, 2;
      $y = 1 / $sp;
      $spe = zeroes($sp) unless defined $spe;
      $ye = $spe / $sp**2;
      PointPlot($w, $x, $y($sel), YERRORBAR=>2*$ye($sel));
      
      #wcols $sel+1, $x, $y($sel);
      
      print "max=", max($y), "\n";
    }
    else
    {
      $y = $h{$cond}{$norm2};
      PointPlot($w, $x, $y, SYMBOLSIZE=>0.9);
    }
    LinePlot($w, pdl(0,10), pdl(0,10), COLOR=>'RED');
    
    
    
    TopText($w, $cond);
  }
  $w->close();
}
