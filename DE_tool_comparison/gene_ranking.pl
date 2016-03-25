#!/sw/bin/perl

=head1 NAME

B<gene_ranking.pl>

=head1 DESCRIPTION

Creates a table of gene ranking according to DE tests. Genes from each DE test are ranked, a score for each gene is calculated and a table (sorted by the score) with all genes from all tools is produced.

There are two ways of ranking genes from individual DE tools (controlled by C<-sig> option). The first method ranks them by the increasing p-values. Where p-values are identical (e.g. equal zero), ranking is (optionally, see option C<-pfc>) done by the decreasing absolute fold change. The score for a gene is then the rank product.

The second method looks only at significance of each gene, according to a given criterion. 0/1 is assigned for significant/non-significant gene. The score is the sum of these 0/1s.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;

use Stats;
use Distribution;
use Tools;
use  GRNASeq;

$| = 1;
Stamp();

my $testinfofile = "de_tests.txt";
my $outfile = 'gene_ranking.csv';
my $tests = 't,lt,mw,perm,bs,bayseq,cuffdiff,degseq,deseq1,deseq2,ebseq,edger,edgerglm,limma,noiseq,poissonseq,samseq';
my $significance;
my $alpha = 0.05;
my $multicor = 'bh';
my $nrep;
my $exclude;
my $allrep;
my $genlist = 'genlist.tsv';
my $pfc;
my ($help, $man);
GetOptions(
  #'nrep=i' => \$nrep,
  #allrep => \$allrep,
  'testinfofile=s' => \$testinfofile,
  'outfile=s' => \$outfile,
  'alpha=f' => \$alpha,
  'multicor=s' => \$multicor,
  #'exclude=s' => \$exclude,
  'genlist=s' => \$genlist,
  'tests=s' => \$tests,
  sig => \$significance,
  pfc => \$pfc,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;


my %ex = map {$_ => 1} split /,/, $exclude if defined $exclude;
my @genes = ReadGeneList($genlist);
my %geneidx = map {$genes[$_] => $_} (0 .. @genes-1);
my %tests = ReadTestfileInfo($testinfofile);
my @tests = split(/,/, $tests);

if($allrep)
{
  for my $n (3 .. 40)
  {
    #next if $n == 10 || $n == 17;
    #print "$n\n";
    my ($ranks, $names) = ReadDEFiles($n, 1);
    my $s = ($significance) ? '_sig' : '';
    my $outf = "./generank/ranking${s}_n$n.csv";
    WriteRanks($outf, $ranks, $names);
  }
  print "\n";
}
else
{
  my ($ranks, $names) = ReadDEFiles($nrep);
  WriteRanks($outfile, $ranks, $names);
}

#################################

sub WriteRanks
{
  my ($file, $ranks, $names) = @_;
  
  my $pvfile = $file;
  $pvfile .= '.pv';
  
  my $tot;
  if($significance)
    {$tot = sumover $ranks}   # sum of 0/1 significances
  else
    {($tot) = exp(statsover(log($ranks)))}   # rank product 
  
  my $idx = qsorti $tot;
  
  my $head = "rank\tgene\tscore\t" . join("\t", @$names);
  my $pvhead = "gene\t" . join("\t", @$names);
  $pvhead =~ s/ /_/g;
  
  local *F;
  local *P;
  open F, ">$file" or die "Cannot open $file for writing\n";
  open P, ">$pvfile" or die "Cannot open $pvfile for writing\n";
  print F "$head\n";
  print P "$pvhead\n";
  my $n = 1;
  for my $k (0 .. @genes - 1)
  {
    my $i = at($idx, $k);
    my $gene = $genes[$i];
    my $r = $ranks(,$i;-);
    #my $tot = sum($r);
    my $p = prodover($r);
    next if at($p, 0) == 0 && !$significance;   # reject zeroes
    
    printf F "%d\t%s\t%8.2f", $n++, $gene, at($tot, $i);
    printf F "\t%7.0f", at($r, $_) for (0 .. nelem($r)-1);
    print F "\n";

    printf P "%s", $gene;
    printf P "\t%7.0f", at($r, $_) for (0 .. nelem($r)-1);
    print P "\n";
  }
}

#################################

sub ReadDEFiles
{
  my ($nrep, $quiet) = @_;
  
  my @files = ();
  my @names = ();
  my @r = ();

  for my $test (@tests)
  {
    next if defined $ex{$test};
    #my ($name, $bsform, $sp_bsform, $fullfile, $cor, $fcsig) = @{$tests{$test}};
    print "$test\n";
    my %t = %{$tests{$test}};
    my $file;
    if($nrep)
      {$file = sprintf $t{repfile}, $nrep}
    else
      {$file = $t{fullfile}}
    print "$file" unless $quiet;
    if($file eq '-')
      {print "  ### skipping non-specified file\n" unless $quiet; next}
    unless(-e $file)
      {print "  ### file $file does not exist!\n"; next}
    #print "\n";
    my $r = ReadDEFile($file, $t{name}, $t{adjustment});
    push @r, $r;
    push @names, $t{name};
  }
  return transpose(pdl(@r)), \@names
}

#################################

sub ReadDEFile
{
  my ($file, $test, $mcor) = @_;
  
  #my ($fc, $p) = rcols $file, 1, 2;
  #my $rank = rank($p);
  
  my @fc = ();
  my @p = ();
  my @gns = ();
  
  local *F;
  open F, $file or die "\nCannot open $file\n";
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    next if $line =~ /^gene/i;
    my ($gene, $fc, $p) = split /\t/, $line;
    $gene =~ tr/A-Z/a-z/;
    $p = 1 if $p =~ /NA/;
    $fc = 0 if $fc =~ /NA/;
    push @p, $p;
    push @fc, $fc;
    push @gns, $gene;
  }
  
  my $fc = pdl(@fc);
  my $p = pdl(@p);
  
  #$fc = log($fc)/log(2) if $file =~ /samseq/;  # silly!
  
  
  
  my $rank;
  
  if($significance)
  {
    #print "MCOR=$mcor\n";
    my $plim = MulticorLimit($p, $mcor, $multicor, $alpha);
    #if($mcor eq 'none')
    #{
    #  if($multicor eq 'hb')
    #    {$plim = HolmBonferroniLimit($p, $alpha)}
    #  elsif($multicor eq 'bh')
    #    {$plim = BenjaminiHochbergLimit($p, $alpha)}
    #  else
    #    {die "Unknown -multicor=$multicor\n"}
    #}
    #elsif($mcor eq $multicor)
    #  {$plim = $alpha}
    #else
    #  {die "Required correction is $multicor, but $test corrected with $mcor\n"}
    $rank = ($mcor =~ /llim/) ? ($p <= $plim): ($p >= $plim);
    
    my $n = nelem($p);
    my $nzer = nelem($p->where($p==0));
    my $nsig = nelem($p->where($p <= $plim));
    
     print "$file\n";
     print "  n = $n  n0 = $nzer  nsig = $nsig  plim = $plim\n";
    #print "$rank\n" if $file =~ /noiseq/;
  }
  else
  {
    $rank = ($pfc) ? rank2($p, $fc) : rank($p);
  }
  
  my $fill = ($significance) ? 1 : 0;
  my @generank = ($fill) x @genes;
  for my $i (0 .. @gns - 1)
  {
    my $gene = $gns[$i];
    my $genei = $geneidx{$gene};
    if(defined $genei)
      {$generank[$genei] = at($rank, $i)}
    else
      {print "Cannot find $gene from $file\n"}
    $i++
  }
  return pdl(@generank)
}


=head1 SYNOPSIS

  gene_ranking.pl -testinfofile=de_tests.dat -outfile=gene_ranking.txt   
    
=head1 OPTIONS

=over 4

=item B<-testinfofile>=I<pathname>

A file with DE tools metadata. See C<plot_fcp.pl> for details.

=item B<-tests>=<string>

Comma-delimited list of test names (as defined in the file specified by -testinfofile) to be included.

=item B<-outfile>=I<filename>

Name of the output file. The default name is C<gene_ranking.csv>.

=item B<-genlist>=I<filename>

A file with a list of gene names (first column) to use. The default name is C<genlist.tsv>.

=item B<-pfc>

If specified, ranking is done on p-values and fold changes (if not specified, ranking is done only on p-values).

=item B<-sig>

If specified, significance will be used to rank genes as opposed to p-value/fold change.

=item B<-multicor>=[none|BH|HB]

Multiple test correction to apply to tests/tools that return uncorrected raw p-values. The default value is C<BH>.


=back
