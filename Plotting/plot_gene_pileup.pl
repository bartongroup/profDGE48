#!/sw/bin/perl -w

=head1 NAME

B<plot_gene_pileup.pl>

=head1 DESCRIPTION

Plots read pileup for a given gene/locus. It shows the mean pileup across replicates, plus/minus one standard deviation. Additionally, a replicate can be highlighted on top of the mean.

Pileup files need to be created using C<build_pileup.pl>. You also need to run C<norms.pl> first in order to created normalization factor files.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use GRNASeq;
use GetOpt;

use PLgraphs;
use Tools;

$| = 1;

Stamp();
ReadDefs();

my ($gene, $locus);
my $margin = 0.5;
my $cs = 0.5;
my ($maxy, $gridlines, $showrep, $all);
GetMyOptions(
  'gene=s' => \$gene,
  'locus=s' => \$locus,
  #'margin=f' => \$margin,
  #'maxy=f' => \$maxy,
  'showrep=s' => \$showrep,
  #'charsize=f' => \$cs,
  #all => \$all,
  #gridlines => \$gridlines,
);
die "Need -gene or -locus\n" unless defined $gene || defined $locus;
die "Provide valid GFF file\n" unless defined $gfffile && -e $gfffile;
die "Provide valid pileup directory\n" unless defined $pileupdir && -d $pileupdir;

$gene =~ tr/A-Z/a-z/ if $gene;
my ($gendat, $genidx) = ReadGFF($gfffile);


my %opt = (
  clean => $clean,
  exclude => $exclude,
  showrep => $showrep,
  margin => $margin,
  norm => $norm,
  cs => $cs,
  maxy => $maxy,
  gridlines => $gridlines
);


my $w = NewWindow($psfile);
  PlotGenePileup($w, 1, \@conds, $gene, $locus, $gendat, $genidx, %opt);
$w->close();


sub PlotGenePileup
{
  my ($w, $pan, $conds, $gene, $locus, $gendat, $genidx, %opt) = @_; 


  my ($chr, $pos1, $pos2);

  if(defined $gene)
  {
    # data for selected gene
    my ($c, $start, $end, $strand, $gene_name, $exons) = GeneStructure($gene, $gendat);
    my $span = $end - $start;
    $chr = $c;
    $pos1 = floor($start - $margin * $span);
    $pos2 = ceil($end + $margin * $span);
  }
  elsif(defined $locus)
  {
    if($locus =~ /^(.+):(\d+)\-(\d+)$/)
    {
      $chr = $1;
      $pos1 = $2;
      $pos2 = $3;
    }
    else
      {die "Unrecognized locus $locus\n"}
  }
  
  # select genes within pos1-pos2 from index
  my ($starts, $ends, $idx, $genes) = @{$genidx->{$chr}};
  my $i1 = vsearch $pos1, $starts;
  my $i2 = vsearch $pos2, $starts;
  $i1-- if $i1 > 0;
  
  # window position
  my $px = ($pan - 1) % $PL_nx + 1;
  my $py = int(($pan - 1) / $PL_nx) + 1;
  my ($xmin, $xmax, $ymin, $ymax) = PanelCoord($px, $py);
  my $struc = 0.1;
  my $rat = 0.2;
    
  # plot gene structure

  my $genwidth = 0.1;
  my $exonwidth = 0.7;
  my $ys1 = $ymin;
  my $ys2 = $ymin + $struc * ($ymax - $ymin);
  $w->xyplot(pdl(-1), pdl(-1),
    VIEWPORT => [$xmin, $xmax, $ys1, $ys2],
    BOX => [$pos1, $pos2, -1, 1],
    %PL_empty
  );
  my $col = [215, 215, 255];
  
  my @grid = ();
  for my $i ($i1 .. $i2)
  {
    my $thisgene = $genes->[at($idx, $i)];
    my ($ch, $st, $en, $strnd, $gname, $exs) = GeneStructure($thisgene, $gendat);
    $gname ||= $thisgene;
    #print "### $ch, $st, $en, $strnd, $gname\n";
    LinePlot($w, pdl($pos1,$pos2), pdl(0,0), COLOR=>$col);
    genebar($w, $st, $en, $strnd, $pos2-$pos1, $genwidth, $col);
    for my $e (@$exs)
    {
      my ($e1, $e2, $str) = @$e;
      push @grid, ($str eq '+') ? [$e1, $e2] : [$e2, $e1];
      genebar($w, $e1, $e2, $str, $pos2-$pos1, $exonwidth, $col)
    }
    my $mid = ($st + $en) / 2;
    $w->text($gname, CHARSIZE=>$cs, COLOR=>'BLACK', TEXTPOSITION=>[$mid, 0, 0, 0, 0.5]);
  }

  my @showrep = ();
  if($opt{showrep})
  {
    my $errmsg = "-showrep requires format 1,5 or 1,- or -,5\n";
    @showrep = split /,/, $opt{showrep};
    die $errmsg unless @showrep == 2;
  }
  
  # plot pileup for conditions

  $PL_deltay = 0.05;
  my $ncond = @conds;
  my $cspan = ($ymax - $rat*($ymax - $ymin) - $ys2 - $PL_deltay) / $ncond;
  my $yc1 = $ys2 + $PL_deltay;
  my $i = 1;
  my @M = ();
  my @S = ();
  my $P;
  my $k = 0;
  for my $cond (@$conds)
  #my $cond = 'WT';
  {
    undef $opt{show};
    $opt{show} = $showrep[$k] if defined $showrep[$k] && !($showrep[$k] eq '-');
    $k++;
    
    my ($pos, $m, $s, $show, $d) = ReadPileup($cond, $chr, $pos1, $pos2, %opt);
    
    my $sel = which(($pos >= $pos1) & ($pos <= $pos2));
    $pos = $pos($sel);
    $m = $m($sel);
    $s = $s($sel);
    $show = $show($sel) if defined $show;
    $d = $d(,$sel);

    $P = $pos;
    push @M, $m;
    push @S, $s / sqrt(nelem($m));
    my $max = max($m + $s) * 1.05;
    $max = $opt{maxy} if defined $opt{maxy};
    $max = max($d) if $all;
    print "Selected ", nelem($pos), " bp\n";
    
    my $cstr = ($chr =~ /ERCC/) ? '' : ' chromosome';
    my $yc2 = $yc1 + $cspan;
    my $xbox = ($i == 1) ? 'BCFNST' : 'BCST'; 
    my $xlab = ($i == 1) ? "Position along$cstr $chr" : '';
    $w->xyplot(pdl(-1), pdl(-1),
      VIEWPORT => [$xmin, $xmax, $yc1, $yc2],
      BOX => [$pos1, $pos2, 0, $max],
      XBOX => $xbox, YBOX => 'BCNST',
      XLAB => $xlab, YLAB => 'Read count',
      MINTICKSIZE => 0.4, MAJTICKSIZE => 0.4,
      CHARSIZE => $cs
    );
    
    if($all)
    {
      my ($ncol, $nrow) = dims($d);
      for my $i (0 .. $ncol - 1)
        {LinePlot($w, $pos, $d($i,;-), COLOR=>'GREY')}
    }
    else
    {
      LinePlot($w, $pos, $m-$s, LINEWIDTH=>1, COLOR=>'GREY');
      LinePlot($w, $pos, $m+$s, LINEWIDTH=>1, COLOR=>'GREY');
    }

    LinePlot($w, $pos, $m, LINEWIDTH=>3);
    LinePlot($w, $pos, $show, LINEWIDTH=>3, COLOR=>'RED') if defined $show;
    gridlines($w, \@grid, 0, $max) if $opt{gridlines};
        
    TopText($w, $cond, x=>0.01);
    TopText($w, $opt{show}, x=>0.99, c=>1, COLOR=>'RED') if $opt{show};
    
    $yc1 = $yc2;
    $i++
  }
  
  return;
  
  # ratio of pileups
  
  my ($m1, $m2) = @M;
  my ($e1, $e2) = @S;
  
  my $zer = which(($m1 == 0) | ($m2 == 0));
  my $sel = which(($m1 > 0) & ($m2 > 0));

  $m1($zer) .= 1;
  $m2($zer) .= 1;
  my $r = log($m2 / $m1) / log(2);
  my $e = sqrt(($e1/$m1)**2 + ($e2/$m2)**2);
  my $max = max(abs($r));
  print "max = $max\n";
  
  $w->xyplot(pdl(-1), pdl(-1),
    VIEWPORT => [$xmin, $xmax, $ymax - ($ymax-$ymin)*$rat, $ymax],
    BOX => [$pos1, $pos2, -$max, $max],
    XBOX => 'BCST', YBOX => 'BCNST',
    XLAB => '', YLAB => 'log2 R',
    MINTICKSIZE => 0.4, MAJTICKSIZE => 0.4,
    CHARSIZE => $cs
  );
  
  LinePlot($w, pdl($pos1, $pos2), pdl(0, 0), COLOR=>'GREY');
  LinePlot($w, $P($sel), $r($sel));
  gridlines($w, \@grid, -$max, $max) if $opt{gridlines};
  
  #LinePlot($w, $P($sel), $r($sel)+$e($sel), COLOR=>'GREY');
  #LinePlot($w, $P($sel), $r($sel)-$e($sel), COLOR=>'GREY');
  
}

sub gridlines
{
  my ($w, $grid, $min, $max) = @_;
  
  for my $g (@$grid)
  {
    my ($g1, $g2) = @$g;
    LinePlot($w, pdl($g1,$g1), pdl($min,$max), COLOR=>[180,255,180]);
    LinePlot($w, pdl($g2,$g2), pdl($min,$max), COLOR=>[255,180,180])
  }
}

sub genebar
{
  my ($w, $p1, $p2, $strand, $totlen, $wd, $col) = @_;
  
  my $len = $p2 - $p1;
  my $tip = $totlen * 0.02;
  $tip = 0.5 * $len if $tip > 0.5 * $len;
  my $pt = $p2 - $tip;

  #BarPlot2($w, pdl($p1), pdl($p2), pdl($width), y0=>-$width, colour=>$col);
  #BarPlot2($w, pdl($p2 - $tip), pdl($p2), pdl($width), y0=>-$width, colour=>[100,100,100]);
  
  my $R = pdl($col->[0], 0);
  my $G = pdl($col->[1], 0);
  my $B = pdl($col->[2], 0);
  plscmap1($R, $G, $B);
  
  my $px = pdl($p1, $pt, $p2, $pt, $p1, $p1);
  my $py = pdl(-$wd, -$wd, 0, $wd, $wd, -$wd);
  
  $px = $p2 - $px + $p1 if $strand eq '-';
  
  plcol1(0);
  plfill($px, $py);
}


sub ReadPileupStats
{
  my ($cond, $chr, $pos1, $pos2) = @_;
  
  my $pfile = "${pileupdir}/${cond}_chr${chr}_pileup_stats.tsv";
  die "Cannot find file $pfile\n" unless -e $pfile;
    
  print "Reading pileup file...";
  my $lines = $pos1-1 . ":" . $pos2;
  my ($pos, $m, $s) = rcols $pfile, {LINES=>$lines};
  print " done\n";
  
  return $pos, $m, $s
}


sub ReadPileup
{
  my ($cond, $chr, $pos1, $pos2, %opt) = @_;
  
  $pos1 = 1 if $pos1 < 1;
  
  my $pfile = "${pileupdir}/${cond}_chr${chr}_pileup.tsv";
  die "Cannot find file $pfile\n" unless -e $pfile;
  
  my $nf = rcols "${cond}_$opt{norm}_normfac.txt", 1 if defined $opt{norm} && $opt{norm} ne 'none';
  
  my ($cl1, $cl2) = split /:/, $CleanExclude;
  my %cl = (WT => $cl1, Snf2 => $cl2);  
  $opt{exclude} = $cl{$cond} if $opt{clean};
  my $ex = pdl(split /,/, $opt{exclude}) - 1 if defined $opt{exclude};
  my $in = pdl(split /,/, $opt{include}) - 1 if defined $opt{include};
  
  local *F;
  open F, $pfile or die;
  my $line = <F>;
  my @s = split " ", $line;
  my $ncol = @s;
  my $lines = $pos1-1 . ":" . $pos2;
  
  print "Reading pileup file $pfile...";
  my ($pos, @d) = rcols $pfile, (0 .. $ncol-1), {LINES=>$lines}; 
  print " done\n";
    
  my $d = pdl(@d)->transpose();
  my ($nrep, $npos) = dims($d);
  print "$nrep x $npos\n"; 
    
  $d /= $nf if defined $nf;

  my $show;
  if(defined $opt{show})
  {
    $show = $d($opt{show}-1,;-)
  }
  
  if(defined $ex)
  {  
    my $all = sequence($nrep);
    my $xsel = setops($all, 'XOR', $ex);
    $d = $d($xsel,);
  }
  if(defined $in)
  {
    $d = $d($in,)
  }
  ($nrep, $npos) = dims($d);
  print "$nrep x $npos\n"; 
  
  
  my ($m, $s) = statsover($d);
  
  return $pos, $m, $s, $show, $d;
}




=head1 SYNOPSIS

  plot_gene_pileup.pl -gene=yhr215w -clean
  plot_gene_pileup.pl -locus=VIII:551800-553500 -clean
      
=head1 OPTIONS

=over 4

=item B<-gene>=I<string>

Gene to plot. If not defined, C<-locus> must be supplied instead.

=item B<-locus>=I<string>

Locus to plot. Format is C<chromosome:pos1-pos2>.

=item B<-clean>

If specified, only clean replicates will be used, as defined in C<defs.dat>

=item B<-norm>=I<string>

Normalization to be used. Needs C<*normfac.txt> files created with script C<norms.pl>. Default value is C<deseq>.

=item B<-showrep>=I<string>

Which replicate to show on top of mean and standard deviation of other replicates. The format is [r1],[r2] for condition 1 and 2. For example 21,6 will highlight replicate 21 in the first condition and 6 in the second one. You can also specify one replicate: 21,- or -,6.

=item B<-psfile>=I<pathname>

Output postscript file.

=back
