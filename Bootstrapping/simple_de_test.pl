#!/sw/bin/perl

=head1 NAME

B<simple_de_test.pl>

=head1 DESCRIPTION

Simple differential expression test. Gene expression data should be stored in two tables, one for each condition. These can be created by the script C<combine_replicates.pl>.

This script can perform a few flavours of t-test, Mann-Whitney test and Kolmogorov-Smirnov test. See option C<-test> for details.

Results are stored in a tab-delimited output file. Columns are:

  * gene name
  * log2 fold change
  * p-value
  * mean for first condition
  * mean for second condition

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
use HTMLTable;
use Stats;

$| = 1;

my $test = 't';
my $outfile;
my $randcond;
my $defsfile;
GetMyOptions(
  'test=s' => \$test,
  'outfile=s' => \$outfile,
  'randcond=s' => \$randcond,
  'defsfile=s' => \$defsfile
);

###############################################

warn "WARNING: KS test does not work properly for discrete data\n" if $test eq 'ks';


ReadDefs($defsfile);
my %data = GetNormalizedData(
  include=>$include,
  nonzero=>$nonzero,
  norm=>$norm,
  type=>$type,
  clean=>$clean,
  spclean=>$spclean,
  randrep=>$randrep,
  randrep2=>$randrep2,
  condition => $randcond
); 

my @myconds = @conds;
if(defined $randrep2)
{
  die "Need -randcond for randrep2\n" unless defined $randcond;
  @myconds = ("$randcond.1", "$randcond.2");
}

my @genes;
my %g2i = ();
my $i = 0;
for my $c (@myconds)
{
  my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel) = @{$data{$c}}; 
  my %g = GenHash($genes, $sel);
  $g2i{$c} = \%g;
  @genes = @$genes if $i == 0;
  $i++ 
}

my ($cond1, $cond2) = @myconds;
my $d1 = $data{$cond1}[0];
my $d2 = $data{$cond2}[0];

my $cl = ($clean) ? '_clean' : '';
my $nm = $norm;
my $rr = ($randrep) ? "_rnd$randrep" : '';
$outfile ||= "${cond1}_${cond2}_${nm}${cl}${rr}_${test}_test.tsv";
open F, ">$outfile";


# prepare data
print "Preparing data...";
my (@dat, @x1, @x2);
for my $gene (@genes)
{
  $gene =~ tr/A-Z/a-z/;
  my $i1 = $g2i{$cond1}{$gene};
  my $i2 = $g2i{$cond2}{$gene};
  if(defined $i1 && defined $i2)
  {
    my $x1 = $d1(,$i1;-);
    my $x2 = $d2(,$i2;-);
    
    #print "$x1\n$x2\n" if $gene eq 'rdn5-1';
    
    my $m1 = mean($x1);
    my $m2 = mean($x2);
    my $fc = log($m2 / $m1) / log(2);
    
    push @dat, [$gene, $m1, $m2, $fc];
    push @x1, $x1;
    push @x2, $x2;
  }
}
print " done\n";
my $n = @dat;


# perform test
print "Performing test...";
my $Ps = ShrinktTest(pdl(@x1), pdl(@x2)) if $test eq 'st';
for my $i (0 .. $n - 1)
{
  my ($gene, $m1, $m2, $fc) = @{$dat[$i]};
  my $x1 = $x1[$i];
  my $x2 = $x2[$i];
 
  my $p;
  if($test eq 't')
    {($p) = tTest($x1, $x2)}
  elsif($test eq 'lt')
    {$p = logtTest($x1, $x2)}
  elsif($test eq 'lnt')
    {$p = logntTest($x1, $x2)}
  elsif($test eq 'mw')
    {$p = MannWhitney($x1, $x2, 1)}
  elsif($test eq 'ks')
    {(my $D, $p) = KStwo($x1, $x2)}
  elsif($test eq 'st')
    {$p = at($Ps, $i)}
  #my ($p) = tTest($x1, $x2);
  printf F "%s\t%.4g\t%.4g\t%.4g\t%.4g\n", $gene, $fc, $p, $m1, $m2;
}

print " done\n";
print "\nResults in $outfile\n";




=head1 SYNOPSIS

  simple_de_test.pl -norm=deseq -clean -test=mw

=head1 OPTIONS

=over 4

=item B<-test>=I<[t|lt|st|mw|ks]>

Type of differential expression test.

  t - standard t-test for two samples with identical variances
  lt - log-ratio t-test using approximately normal statistic Z = (log m2 - log m1) / sqrt(((s1/m1)**2)/n1 + ((s2/m2)**2)/n2), where m1, m2 are sample means and s1, s2 are sample standard deviations
  st - shrinkage variance test of Opgen-Rhein and Strimmer (2007)
  mw - Mann-Whitney test
  ks - Kolmogorov-Smirnov test; B<warning> - KS test does not work properly for discrete data!
  
=item B<-outfile>=I<string>

Output file. If not defined, a default name will be build using condition names and type of test.

=item B<-norm>=I<string>

Expression data normalization. Input expression files contain raw counts. These are normalized on the fly, using one of the following methods:

  none - no normalization
  deseq - DESeq normalization (default)
  tmm - trimmed mean of M values normalization
  fpkm - approximate FPKM normalization (not identical to cuffdiff!)
  totcount - total count
  totcountm - total count with normalization factors scaled to the mean of 1

=item B<-clean>

A switch indicating that only clean replicates should be used. Clean replicates are defined in C<GRNASeq.pm> module.

=item B<-defsfile>=I<string>

Path to defs.dat configuration file.

=item B<-randcond>=I<string>

If specified, test will be performed on randomly selected replicates from the same condition, defined here.

=back
