#!/sw/bin/perl

=head1 NAME

B<gene_ranking.pl>

=head1 DESCRIPTION

Creates a table of gene ranking according to DE tests. Genes from each DE test are ranked, a score for each gene is calculated and a table (sorted by the score) with all genes from all tools is produced.

This script uses a different approach to C<gene_ranking.pl>.  It looks at a proportion of significant DE calls across bootstraps for each gene. This is calculated by C<make_powerstats_db.pl> and stored in files *_sigprop.stats. This scripts simply reads these files and collates information.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;

use Tools;
use Stats;
use GRNASeq;

$| = 1;
Stamp();

my $testinfofile = "de_tests.txt";
my $dir = 'powerstats_db_ref';
my $form = 'power_%s_true-%s_sigprop.stats';
my $tests = 'lt,edger,limma,deseq,cuffdiff,bayseq,samseq,noiseq,degseq,poissonseq';
my $nrep;
my $outfile = 'gene_ranking_sigprop.csv';
my $genlist = 'genlist.tsv';
my ($help, $man);
GetOptions(
  'dir=s' => \$dir,
  #'form=s' => \$form,
  'nrep=i' => \$nrep,
  'tests=s' => \$tests,
  'outfile=s' => \$outfile,
  'genlist=s' => \$genlist,
  'testinfofile=s' => \$testinfofile,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;


my @genes = ReadGeneList($genlist);
my $ngen = @genes;
my %geneidx = map {$genes[$_] => $_} (0 .. @genes-1);
my @tests = split(/,/, $tests);
my %tests = ReadTestfileInfo($testinfofile);

if(defined $nrep)
{
  my ($sigprop, $tst) = ReadSigpropFiles($nrep);
  my ($ncol, $nrow) = dims $sigprop;
  print "$ncol x $nrow :", scalar @$tst, "\n";
  my @names = map {$tests{$_}{name}} @$tst;
  WriteRanks($outfile, $sigprop, \@names);
}
else
{
  for my $n (2 .. 40)
  {
    #next if $n == 10 || $n == 17;
    print "$n ";
    my ($sigprop, $tst) = ReadSigpropFiles($n, 1);
    my $outf = "./generank_sigprop/ranking_n$n.csv";
    my @names = map {$tests{$_}{name}} @$tst;
    WriteRanks($outf, $sigprop, \@names);
  }
  print "\n";
}

#################################

sub WriteRanks
{
  my ($file, $sigprop, $names) = @_;
  
  my ($tot) = statsover($sigprop);
  #print "NT = ", nelem($tot), "\n";
  my $idx = qsorti $tot;
  $idx = $idx(-1:0);
  
  my $head = "rank\tgene\tscore\t" . join("\t", @$names);
  
  local *F;
  open F, ">$file" or die;
  print F "$head\n";
  my $n = 1;
  for my $k (0 .. @genes - 1)
  {
    my $i = at($idx, $k);
    my $gene = $genes[$i];
    my $r = $sigprop(,$i;-);
    
    printf F "%d\t%s\t%.4f", $n++, $gene, at($tot, $i);
    printf F "\t%.2f", at($r, $_) for (0 .. nelem($r)-1);
    print F "\n";
  }
}

#################################

sub ReadSigpropFiles
{
  my ($nrep, $quiet) = @_;
  
  local *F;
  
  my @tst = ();
  my @dat = ();
  for my $test (@tests)
  {
    my $file = sprintf "%s/$form", $dir, $test, $test;
    print "$file" unless $quiet;
    unless(-e $file)
      {print "  ### file $file does not exist!\n"; next}
    
    my $sigprop = zeroes($ngen);  
    open F, $file or die;
    <F>;
    while(my $line = <F>)
    {
      chomp $line;
      my ($n, $gene, $fc, $sp) = split /\t/, $line;
      next unless $n == $nrep;
      $sigprop($geneidx{$gene}) .= $sp
    }
    close F;
    push @dat, $sigprop;
    push @tst, $test; 
    print "\n" unless $quiet;
  }
  my $D = pdl(@dat)->transpose();
  
  return $D, \@tst
}


=head1 SYNOPSIS

  gene_ranking.pl -testinfofile=de_tests.dat -outfile=gene_ranking.txt   
    
=head1 OPTIONS

=over 4

=item B<-dir>=I<path>

Directory with results from C<make_powerstats_db.pl>. 

=item B<-testinfofile>=I<pathname>

A file with DE tools metadata. See C<plot_fcp.pl> for details.

=item B<-outfile>=I<filename>

Name of the output file. The default name is C<gene_ranking.csv>.

=item B<-genlist>=I<filename>

A file with a list of gene names (first column) to use. The default name is C<genlist.tsv>.

=item B<-nrep>=I<number>

Use the specific number of replicates (from bootstrap tests) to do ranking. If not specified, the full clean replicate set results (as specified in defs.dat file) will be used.

=item B<-allrep>

Loop through all replicates and produce one ranking file for each.

=item B<-pfc>

If specified, ranking is done on p-values and fold changes (if not specified, ranking is done only on p-values).

=item B<-sig>

If specified, significance will be used to rank genes as opposed to p-value/fold change.

=item B<-multicor>=[none|BH|HB]

Multiple test correction to apply to tests/tools that return uncorrected raw p-values. The default value is C<BH>.

=item B<-exclude>=I<list>

A comma delimited list of tools/tests to exclude from analysis.

=back
