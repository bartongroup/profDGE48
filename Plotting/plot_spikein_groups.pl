#!/sw/bin/perl

=head1 NAME

B<>

=head1 DESCRIPTION


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

use CompBio::PLgraphs;
use CompBio::Tools;
use CompBio::Distribution;

$| = 1;
Stamp();

my $file;
my $maxrep = 26;
my $fmax = 3;
my $testinfofile = "de_tests.txt";
my $genlist = 'goodspikes.txt';
my $sellist = 'spikes_unregulated.txt,spikes_down50.txt,spikes_down67.txt,spikes_up.txt';
my $fclines = '0.5,0.67,1,4';
my $psfile;
GetOptions(
  'testinfofile=s' => \$testinfofile,
  'genlist=s' => \$genlist,
  'sellist=s' => \$sellist,
  'fclines=s' => \$fclines,
  'psfile=s' => \$psfile,
  'maxrep=i' => \$maxrep
);


###############################################

my @fclines = ();
@fclines = split ',', $fclines if defined $fclines;

die "Need -testinfofile\n" unless defined $testinfofile && -e $testinfofile;
my %alltests = ReadTestfileInfo($testinfofile);
my @tests = ();
for my $test (sort keys %alltests)
{
    my %t = %{$alltests{$test}};
    next if $t{sp_repfile} eq '-';
    my $spfile = sprintf $t{sp_repfile}, $maxrep;
    push @tests, $test if -e $spfile;
}

my $n = @tests;
my $w = NewWindow($psfile);

my $div = 0.7;
$PL_ymin = 0.05;
$PL_ymax = $div;

PlotPanelN($w, 1,
  BOX => [0, $n+1, -$fmax, $fmax],
  XBOX => '',
  YBOX => 'BNST',
  YLAB => 'log#d2#u FC'
);

LinePlot($w, pdl(0, $n+1), pdl(log($_)/log(2), log($_)/log(2)), COLOR=>'GREY') for (@fclines);

my $pos = 1;
for my $test (@tests)
{
  my %t = %{$alltests{$test}};
  my $spfile = sprintf $t{sp_repfile}, $maxrep;
  print "$spfile ";
  my ($dat, $M, $S) = GetMeanSpikes($spfile, $t{FCsig}, $genlist);
  print "$M $S\n";
  my $y = zeroes($M) + $pos;
  #PointPlot($w, $M, $y, XERRORBAR=>2*$S);
  PlotPanelN($w, 1, BOX=>[-$pos, -$pos+$n+1, -$fmax, $fmax], %PL_empty);
  BoxPlot($w, [$_], boxwidth=>0.2) for (@$dat);
  $pos++;
} 

$PL_ymin = $div;
$PL_ymax = 0.95;

PlotPanelN($w, 1,
  BOX => [0, $n+1, 0, 1],
  %PL_empty
);


$pos = 1;
for my $test (@tests)
{
  $w->text($alltests{$test}{name}, CHARSIZE=>0.6, COLOR=>'BLACK', TEXTPOSITION=>[$pos++, 0, 0, 1, 0]);
} 

$w->close();


###########################################

sub GetMeanSpikes
{
  my ($file, $fcsign, $genlist) = @_;
  
  my ($p, $lfc, $genes) = ReadFCP($file, $fcsign, $genlist);
  my %g2i = GenHash($genes);
  #$lfc *= $fcsign;
  
  my @d = ();
  my @M = ();
  my @S = ();
  my @selfiles = split /\,/, $sellist;
  for my $selfile (@selfiles)
  {
    die "Cannot find file $selfile\n" unless -e $selfile;
    my @selgen = ReadGeneList($selfile);
    my @s = map {$g2i{$_}} @selgen;
    my $sel = pdl(@s);
    my $dat = $lfc($sel);
    my ($mean, $sd) = stats($dat);
    #print "$selfile M=$mean S=$sd\n";
    #print "$sel   ", $lfc($sel), "\n";
    push @M, $mean;
    push @S, $sd;
    push @d, $dat;
  }

  return \@d, pdl(@M), pdl(@S);
}

=head1 SYNOPSIS

  plot_fcp.pl -file=de_test_results.csv     
    
=head1 OPTIONS

=over 4


=item B<-file>=I<filename>

A file with DE test results. As a minimum, it should contain a column with p-values and a column with log2 fold change. By default, the first column is a gene name, the second column is log2 fold change and the third column is the p-value. These can be changed by column options described below.

=item B<-colg>=I<number>

Column with gene name (identifier). The default value is 1.

=item B<-colfc>=I<number>

Column with log2 fold change. The default value is 2.

=item B<-colp>=I<number>

Column with p-value. The default value is 3.

=item B<-fcsign>=[-1|1]

Sign of log2 fold change. If the tool reports fold change the worng way around, use -1 to reverse it.

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

=item B<-sellist>=I<filename>

An optional comma-delimited list of files, containing gene names to highlight (in the first column). Genes from each file will be highlighed in a different colour.

=item B<-genlist>=I<filename>

An optional file with gene names (first column) to be displayed. Other genes will be ignored.

=back
