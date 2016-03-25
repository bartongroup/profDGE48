#!/sw/bin/perl

=head1 NAME

B<plot_fcp.pl>

=head1 DESCRIPTION

Plot p-value versus fold-change. The input file for this script is a result from a DE tool, which should be a tab-delimited file containing a column with p-values and a column with log2 fold-change. A column with gene names allows for gene selection and highlighting.

=cut


use strict;
use warnings;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use Math::CDF qw(:all);

use GRNASeq;
use GetOpt;

use PLgraphs;
use Tools;
use Distribution;

$| = 1;
Stamp();

my $file;
my $alpha = 0.05;
my ($bigger, $corrected);
my $title = 'DE test results';
my $lmin;
my $fmax = 3;
my $colg = 1;
my $colfc = 2;
my $colp = 3;
my $fcsign = 1;
my $all;
my $testinfofile = "de_tests.txt";
my $sellist;
my $fclines;
my $label;
my $mcor = 'none';
GetMyOptions(
  'title=s' => \$title,
  'file|f=s' => \$file,
  'alpha=f' => \$alpha,
  'lmin=f' => \$lmin,
  'fmax=f' => \$fmax,
  'colg=i' => \$colg,
  'colfc=i' => \$colfc,
  'colp=i' => \$colp,
  'fcsign=i' => \$fcsign,
  'testinfofile=s' => \$testinfofile,
  'sellist=s' => \$sellist,
  'fclines=s' => \$fclines,
  'label=s' => \$label,
  'mcor=s' => \$mcor,
  all => \$all
);
die "Need -file or -all\n" unless defined $file || $all;


###############################################

my @fclines = ();
@fclines = split ',', $fclines if defined $fclines;

if($all)
{
  die "Need -testinfofile\n" unless defined $testinfofile && -e $testinfofile;
  my %tests = ReadTestfileInfo($testinfofile);
  for my $test (keys %tests)
  {
    print "Plotting $test...";
    #my ($name, $bsform, $sp_bsform, $fullfile, $cor, $fcsig) = @{$tests{$test}};
    my %t = %{$tests{$test}};
    
    my $tit = "DE test results for $t{name}";
    my $ps = "./fcp/fcp_${test}.ps";
    undef $lmin;
    PlotFCP($t{fullfile}, $tit, $t{adjustment}, $t{FCsig}, $ps);
    print "\n";
  }
}
else
  {PlotFCP($file, $title, $mcor, $fcsign, $psfile, $genlist, $sellist, \@fclines)}


###############################################

sub PlotFCP
{
  my ($file, $title, $mcor, $fcsign, $psfile, $genlist, $sellist, $fclines) = @_;
  
  my ($p, $lfc, $genes) = ReadFCP($file, $fcsign, $genlist);
  my %g2i = GenHash($genes);

  
  my @hl = ();                       # highlight
  if(defined $sellist)
  {
    my @selfiles = split (/\,/, $sellist);
    for my $selfile (@selfiles)
    {
      die "Cannot find file $selfile\n" unless -e $selfile;
      my @selgen = ReadGeneList($selfile);
      my @s = map {$g2i{$_}} @selgen;
      push @hl, pdl(@s);
    }
  }
  
  my $n = nelem($p);
  my $n_plus = nelem($lfc->where($lfc>=0));
  my $n_minus = nelem($lfc->where($lfc<0));
  #print "$n $n_plus $n_minus\n";
  
   
  my $lim = MulticorLimit($p, $mcor, $multicor, $alpha);
  print "Limit = $lim\n";
  my $ns = nelem($p->where($p<$lim));
  my $ns_plus = nelem($p->where(($p<$lim) & ($lfc>=0)));
  my $ns_minus = nelem($p->where(($p<$lim) & ($lfc<0)));
  
  print "N = $n, Ns = $ns, N_up = $n_plus, N_down = $n_minus, Ns_up = $ns_plus, Ns_down = $ns_minus\n";
  
  
  my $min = (defined $lmin) ? 10**$lmin : min($p->where($p>0));
  my $zero = $min * 0.5;
  $p->where($p <= $min) .= $zero;
  $lim = ($lim > 0) ? log10($lim) : log10($zero);
  #print "min=$min\n";
  my $nz = which($p > $min);
  my $z = which($p <= $min);
  
  my $lp = log10($p);
  
  # mark fold change outside the plot
  my ($minx, $maxx) = (-$fmax, $fmax);
  my $fcmax = $maxx - 0.1;
  my $foutp = which($lfc > $fcmax);
  my $foutm = which($lfc < -$fcmax);
  my $fout = append($foutp, $foutm);
  my $fin = which(abs($lfc) <= $fcmax);
  $lfc($foutp) .= $fcmax;
  $lfc($foutm) .= -$fcmax;
  
  my $points_out = append($fout, $z);
  my $points_in = setops(sequence($p), 'XOR', $points_out); 
  
  print "min = $min, zero = $zero\n";
  
  #if(defined $lmin)
   # {$zero = 10**($lmin + 0.4)}
  #else
  #  {$lmin = log10($zero) - 0.4}
  my $nbin = 100;
  my ($x1, $h1) = BuildHistogram($lfc($fin),$nbin,$minx,$maxx,1,1);
  my ($x2, $h2) = BuildHistogram($lp($nz),$nbin,$lmin,0,1,1);
  
  my $w = NewWindow($psfile);
  
  #my ($minx, $maxx) = BoxMinMax($lfc);
  my ($miny, $maxy) = BoxMinMax($lp);
  #$miny = $lmin;
  $miny = log10($zero) - 0.4;
  $maxy = 1e-3;
  
   
  my $div = 0.85;
  my $cmax = 0.93;
  
  ###
  
  my $Z = 1e-3;
  
  FullBox($w);
  TopText($w, $title, x=>0.5, y=>-1, c=>0.5, CHARSIZE => 0.7);
  TopText($w, "Data from $file", x=>0.5, c=>0.5, y=>-3.5, CHARSIZE=>0.5);
  TopText($w, "n = $n, n#ds#u = $ns, n#d+#u = $n_plus, n#d-#u = $n_minus", x=>0.5, c=>0.5, y=>-5, CHARSIZE=>0.6);
  
  ###
  
 
  $PL_xmin = 0.1;
  $PL_xmax = $div;
  $PL_ymin = $div;
  $PL_ymax = $cmax;
  my $max1 = max($h1)*1.02;
  PlotPanelN($w, 1, BOX => [$minx, $maxx, 0, $max1], %PL_empty);
  BarPlot($w, $x1, $h1, colour=>[215,215,255], boxes=>1);
  LinePlot($w, pdl(0,0), pdl(0,$max1), LINEWIDTH=>3);
  PlotPanelN($w, 1,
    BOX => [$minx, $maxx, 0, max($h1)*1.02],
    XBOX => 'BSTI', YBOX => '',
    non => 1
  );
  
  $PL_ymin = 0.1;
  $PL_ymax = $div;
  $PL_xmin = $div;
  $PL_xmax = $cmax;
  PlotPanelN($w, 1, BOX => [0, max($h2)*1.02, $miny, $maxy], %PL_empty);
  BarPlot($w, $x2, $h2, colour=>[215,215,255], boxes=>1, swap=>1);
  PlotPanelN($w, 1,
    BOX => [0, max($h2)*1.02, $miny, $maxy],
    XBOX => '', YBOX => 'BSTI',
    non => 1
  );
  
  $PL_xmin = 0.1;
  $PL_ymin = 0.1;
  $PL_ymax = $div;
  $PL_xmax = $div;
  
  PlotPanelN($w, 1,
    BOX => [$minx, $maxx, $miny, $maxy],
    XLAB => "log#d2#u FC",
    YLAB => "log p",
    CHARSIZE => 0.8
  );
  
  #LinePlot($w, pdl(log($_)/log(2), log($_)/log(2)), pdl($miny, $maxy), COLOR=>'GREY') for (@$fclines);
  LinePlot($w, pdl(log($fclines[$_])/log(2), log($fclines[$_])/log(2)), pdl($miny, $maxy), COLOR=>$colours[$_]) for (0 .. @fclines-1);
  
  my $col = (@hl > 0) ? 'GREY' : 'BLACK';
  
  LinePlot($w, pdl($minx,$maxx), pdl(0,0), COLOR=>'GREY');
  LinePlot($w, pdl(0,0), pdl($miny,0), COLOR=>'GREY') unless $fclines;
  PointPlot($w, $lfc($points_in), $lp($points_in), SYMBOLSIZE=>0.6, COLOR=>$col);
  PointPlot($w, $lfc($points_out), $lp($points_out), SYMBOLSIZE=>0.6, COLOR=>'GOLD2');
  LinePlot($w, pdl(-10,10), pdl($lim,$lim), COLOR=>'RED', LINESTYLE=>4);

  for my $i (0 .. @hl - 1)
  {
    PointPlot($w, $lfc($hl[$i]), $lp($hl[$i]), SYMBOLSIZE=>1.3, COLOR=>$colours[$i]);
    #print $hl[$i], " ", $lfc($hl[$i]), "\n"
  }
  TopText($w, $label, CHARSIZE=>1.2, x=>0.02) if $label;
  
  GetPos($lfc, $lp, $genes) unless $psfile;
  
  $w->close();
}





sub GetPos
{
  my ($xx, $yy, $g) = @_;
  
  my %gin;
  do {
    %gin = plGetCursor();
    if($gin{button})
    {
      my $x = $gin{wX};
      my $y = $gin{wY};
      
      my $d = ($x - $xx)**2 + ($y - $yy)**2;
      my $mini = minimum_ind($d);
      my $gene = $g->[$mini];
      
      my $p = 10**$y;
      my $fc = 2**$x;
      printf "%s: fc = %.3g  p = %.3g\n", $gene, $fc, $p;
    }
  } until 0;
}




=head1 SYNOPSIS

  plot_fcp.pl -file=de_edger.txt -psfile=de_edger.ps     
    
=head1 OPTIONS

=over 4

=item B<-file>=I<filename>

A file with DE test results. As a minimum, it should contain a column with p-values and a column with log2 fold change. By default, the first column is a gene name, the second column is log2 fold change and the third column is the p-value. These can be changed by column options described below.

=item B<-psfile>=I<pathname>

Name of the output postscript file.

=item B<-colg>=I<number>

Column with gene name (identifier). The default value is 1.

=item B<-colfc>=I<number>

Column with log2 fold change. The default value is 2.

=item B<-colp>=I<number>

Column with p-value. The default value is 3.

=item B<-fcsign>=[-1|1]

Sign of log2 fold change. If the tool reports fold change the wrong way around, use -1 to reverse it.

=item B<-mcor>=C<none|bh|hb>

Multiple test correction to apply in the plot.

  none - the default value under the assumption that the tool/test output is already corrected
  bh - Benjamini-Hochberg
  hb - Holm-Bonferroni


=item B<-all>

If used, fold-chage/p-value plots will be created for all tools in the file specified by C<-testinfofile>.

=item B<-testinfofile>=I<pathname>

A file with DE tools metadata. This is a text file contating records for each test/tool. Each record has a following format:

  TEST deseq
  name = DEseq
  repfile = /cluster/gjb_lab/cdr/GRNAseq/analysis/bootstraps/DEseq/deseq_%dreps_defaultnorm.tsv
  sp_repfile = /cluster/gjb_lab/cdr/GRNAseq/analysis/bootstraps_with_spike_ins/deseq/deseq_%dreps_package_default.tsv
  fullfile = wrap_deseq.txt
  adjustment = bh
  FCsig = -1

C<name> is the name of the test to be displayed. C<repflie> is a format for bootstrap results (an FCP file) for all replicates (%d stands for replicate number). C<sp_repfile> is for spike-in bootstrap results. C<fullfile> is the result for full clean replicate set. C<adjustment> is the multiple test adjustment used to obtain p-values. C<FCsig> is an optional sign correction for fold-change.

FCP files are tab-delimited with (at least) three columns: gene id, log2 fold change and p-value.

=item B<-sellist>=I<string>

An optional comma-delimited list of files, containing gene names to highlight (in the first column). Genes from each file will be highlighed in a different colour.

=item B<-genlist>=I<filename>

An optional file with gene names (first column) to be displayed. Other genes will be ignored.

=item B<-fclines>=I<string>

Comma-delimited list of fold-change values (not logarithms!) to be shown as vertical lines.

=item B<-title>=I<string>

Title of the plot.

=item B<-alpha>=I<value>

Significance level. The default value is 0.05. A dashed red line at the significance level will be shown in the plot and the number of the significantly DE genes, n_s, reported.

=item B<-lmin>=I<value>

Minimum value (as a logarithm) for the vertical axis. If not specified, the script will calculate it from the data.


=item B<-fmax>=I<value>

Minimum/maximum value (as a base-2 logarithm) for the horizontal axis. The default value is 3.


=item B<-label>=I<string>

Optional label inside the plot (left-top corner).




=back
