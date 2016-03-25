package GRNASeq;
require Exporter;
require PDL;
require PDL::NiceSlice;
require PDL::Graphics::PLplot;
require Carp;
require Tools;
require PLgraphs;
require Distribution;
require Normality;
require Math::CDF;

=head1 NAME

GRNASeq 

=head1 SYNOPSIS

  use PDL;
  use GRNASeq;
  
=head1 DESCRIPTION
  
=head1 FUNCTIONS

=over 4

=cut

require Exporter;

use strict;
use warnings;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;
use Carp;

use Tools;
use PLgraphs;
use Distribution;
use Normality;

use Math::CDF qw(:all);

our $VERSION = 1.1;

our @ISA = ("Exporter");
our @EXPORT = qw(
  $outlierdir
  $countsdir
  $pileupdir
  ReadDefs
  ReadExpressionFile
  ReadGeneList
  ReadGeneCross
  ReadFCP
  NormalizeData
  CountFile
  OutliersFile
  GetNormalizedData
  PlotDistWithFits
  RejectOutliersSigma
  RejectOutliersChauvenet
  ReadGeneLengths
  ReadGFF
  GeneStructure
  GenHash
  HashList
  CreateDiscreteDist
  NormalityFromData
  PoissonFromData
  NBHGFromData
  NBKSFromData
  ReadDETestFile
  ReadTwoDETestFiles
  ReadSpikeAndDETestFiles
  ReadGOFile
  SpikeinFile
  ReadTestfileInfo
  MulticorLimit
  rank2
  PowerTest
  SpikeinPowerTest
  MultiCorLimit
  $CleanExclude
  $SpikeinCleanExclude
  $gfffile
  @conds
  @colours
);

#our @colours = ('BLUE', 'RED', [50,150,50], [234,193,23], [213,71,153], 'GREY');
our @colours = ([0,0,255], [240,50,50], [50,150,50], [187,154,18], [213,71,153], [190,190,190]);

my $defsfile = 'defs.dat';

our ($topdir, $countsdir, $outlierdir, $pileupdir, $bootstrapdir);
our $gfffile;
our ($CleanExclude, $SpikeinCleanExclude); 

our %defs;
our @conds;

my @required = qw(topdir countsdir outlierdir pileupdir bootstrapdir gfffile clean spclean conditions);

################################################################

sub ReadDefs
{
  my ($indefsfile) = @_;
  local *F;
  if (!defined($indefsfile)) {
  	open F, $defsfile or die "Cannot open $defsfile\n";
  } else {
	open F, $indefsfile or die "Cannot open $indefsfile\n";
  }
  
  %defs = ();
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    if($line =~ /^\s*(\S+)\s*=\s*(\S+)\s*$/)
      {$defs{$1} = $2;}
  }
  
  for my $r (@required)
    {die "$r not defined in $defsfile\n" unless defined $defs{$r}}
  
  @conds = split(/,/, $defs{conditions});
  
  $topdir = $defs{topdir};
  $countsdir = "$topdir/$defs{countsdir}";
  $outlierdir = "$topdir/$defs{outlierdir}";
  $pileupdir = "$topdir/$defs{pileupdir}";
  
  $CleanExclude = $defs{clean};
  $SpikeinCleanExclude = $defs{spclean};
  
  $gfffile = $defs{gfffile};
}

################################################################


=item ReadExpressionFile

  my ($d, $m, $s, $sel) = ReadExpressionFile($file, %options);

  Read data from expression file $file. Expression file should be tab-delimited with the first column containing gene/transcript name, and the remaining columns containing replicates. This script returs a two-dimensional piddle $d (genes in rows and replicates in columns) and one-dimensional piddles $m and $s, containing mean and standard deviation, respectively. The last value is a selection piddle.
 
Options:
  
  nonzero=>n selects only rows (genes) with at least n non-zero replicates.
  
  exclude=>[i1, i2, ...,in] excludes a list of replicates from data.

  ranrep=> selects n replicates at random

=cut


sub ReadExpressionFile
{
  my ($file, %opt) = @_;
  
  local *F;
  open F, $file or croak "\nCannot open $file\n";
  <F>;
  my $test = <F>;
  my @a = split " ", $test;
  my $n = @a;
  
  my @d = rcols $file, (1 .. $n-1);
  my $d = pdl(@d)->transpose();
  my ($nrep, $ngen) = dims($d); 

  my ($sel, $xsel);
  if(defined $opt{exclude})
  {
    my $ex = pdl($opt{exclude}) - 1;
    #print "EX = $ex\n";
    my $all = sequence($nrep);
    $xsel = setops($all, 'XOR', $ex);
    #print "ex=$ex\nall=$all\nsel=$sel\n";die;
    $d = $d($xsel,);
    ($nrep, $ngen) = dims($d)
  }
  elsif(defined $opt{include})
  {
    $xsel = pdl($opt{include}) - 1;
    #print "IN = ", qsort($xsel), "\n";
    $d = $d($xsel,);
    ($nrep, $ngen) = dims($d)
  }
  else
    {$xsel = sequence($nrep)}
    
  if(defined $opt{randrep})
  {
    if($opt{randrep} < $nrep)
    {
      my $r = qsorti random $nrep;
      $r = $r(0:$opt{randrep} - 1);  # random selection
      $xsel = $xsel($r);
      $d = $d($r,);
      #print "### $xsel\n";
    }
  }
  
  if(defined $opt{nonzero})
  {
    $opt{nonzero} = $nrep if $opt{nonzero} < 0;
    my $N = sumover($d>0);
    $sel = which($N >= $opt{nonzero});
    $d = $d(,$sel);
  }
  my ($m, $s) = statsover($d);
  
  return $d, $m, $s, $sel, $xsel;
}

################################################################

=item ReadGenelist

  my @g = ReadGeneList($file);
  
  Reads a list of genes (first column) from expression file.

=cut

sub ReadGeneList
{
  my ($file, $skip) = @_;
  
  $skip ||= 0 ;
  
  local *F;
  open F, $file or die "Cannot open $file\n";
  #print "Reading $file...";
  <F> for (1 .. $skip);
  my @g = ();
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    my @s = split " ", $line;
    my $gene = shift @s;
    #next if @s == 0;
    $gene =~ tr/A-Z/a-z/;
    push @g, $gene
  }
  #print " done\n";
  return @g
}

################################################################

sub ReadTotCountFile
{
  my ($file, %opt) = @_;
  
  return undef unless -e $file;
  my @cnt = ();
  local *F;
  open F, $file or die "Cannot open total count file $file\n";
  while(my $line = <F>)
  {
    chomp $line;
    my ($file, $rep, $count) = split " ", $line;
    push @cnt, $count;
  }
  my $cnt = pdl(@cnt);
  my $nrep = nelem($cnt);

  if(defined $opt{repsel})
  {
    $cnt = $cnt($opt{repsel});
  }
  elsif(defined $opt{exclude})
  {
    my $ex = pdl($opt{exclude}) - 1;
    my $all = sequence($nrep);
    my $sel = setops($all, 'XOR', $ex);
    #print "ex=$ex\nall=$all\nsel=$sel\n";die;
    $cnt = $cnt($sel);
  }
  elsif(defined $opt{include})
  {
    my $sel = pdl($opt{include}) - 1;
    $cnt = $cnt($sel);
  }
  return $cnt;  
}

################################################################

sub ReadSpikeinFile
{
  my ($file, %opt) = @_;
  
  return undef unless -e $file;
  my $r = rcols $file, 1;
  my $f = 1 / $r;
  my $nrep = nelem($f);
  
  
  if(defined $opt{repsel})
  {
    $f = $f($opt{repsel});
  }
  elsif(defined $opt{exclude})
  {
    my $ex = pdl($opt{exclude}) - 1;
    my $all = sequence($nrep);
    my $sel = setops($all, 'XOR', $ex);
    #print "ex=$ex\nall=$all\nsel=$sel\n";die;
    $f = $f($sel);
  }
  elsif(defined $opt{include})
  {
    my $sel = pdl($opt{include}) - 1;
    $f = $f($sel);
  }
  return $f;  
}

#############################################################

=item ReadTwoDETestFiles

  my ($genes, $p1, $p2, $fc1, $fc2) = ReadTwoDETestFiles($file1, $file2);
  
Reads two files with results of DE analysis. Each of these files should contain gene names in the first column, log2 fold changes in the second column and p-values in the third column.

Returned are: array of genes $genes, and two piddles with corresponding p-values, $p1 and $p2 and fold changes $fc1, $fc2. 

=cut

sub ReadTwoDETestFiles
{
  my ($f1, $f2) = @_;
  
  my %h1 = getgenehash($f1);
  my %h2 = getgenehash($f2);
  
  my @p1 = ();
  my @p2 = ();
  my @fc1 = ();
  my @fc2 = ();
  my @g = ();
  for my $g (keys %h1)
  {
    if(defined $h2{$g})
    {
      push @g, $g;
      push @p1, $h1{$g}[1];
      push @p2, $h2{$g}[1];
      push @fc1, $h1{$g}[0];
      push @fc2, $h2{$g}[0];
    }
  }
  return \@g, pdl(@p1), pdl(@p2), pdl(@fc1), pdl(@fc2);
}

=item ReadDETestFile

  my %g = ReadDETestFile($file);
  
Reads a file with results of DE analysis. It should contain gene names in the first column, log2 fold changes in the second column and p-values in the third column.

Returned is a hash: gene => [fc, p] 

=cut

sub ReadDETestFile
{
  getgenehash(@_)
}

sub getgenehash
{
  my ($file) = @_;
  
  local *F;
  my %h;
  open F, $file or die "Cannot open file $file\n";
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    next if $line =~ /^gene/i;
    my ($gene, $fc, $p) = split /\t/, $line;
    
    die "$file: $gene undefined p\n" unless defined $p;
    
    $p = 1 if $p =~ /NA/;
    $fc = 0 if $fc =~ /NA/;
    $gene =~ s/\"//g;
    $gene =~ tr/A-Z/a-z/;
    #print "### $gene\n" unless defined $p;
    $p = 0 if $p eq 'NA';
    $fc = 0 if $fc eq 'NA';
    $h{$gene} = [$fc, $p]
  }
  return %h; 
}


#############################################################

=item ReadSpikeAndDETestFiles

  my ($genes, $sp, $p) = ReadSpikeAndDETestFiles($DEfile, $spikefile)

Reads a DE file (fcp format) and a spike-in file (a list of good spikes with columns: spikein @col1 and expected fold change @col5) and combines them together. Returns a list of gene (spike-in) names, spike-in significance (0 for not changing, 1 for changing) and p-values from the DE file. 

=cut

sub ReadSpikeAndDETestFiles
{
  my ($f, $fspike) = @_;
  
  my %h = getgenehash($f);
  my %hspike = getspikehash($fspike);
  
  my @g = ();
  my @sp = ();
  my @p = ();
  for my $spike (sort keys %hspike)
  {
    if(defined $h{$spike})
    {
      push @g, $spike;
      push @sp, $hspike{$spike};    # significance
      push @p, $h{$spike}[1];
    }
  }
  return \@g, pdl(@sp), pdl(@p) 
}

sub getspikehash
{
  my ($file) = @_;
  
  local *F;
  my %h;
  open F, $file or die "Cannot open file $file\n";
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    next if $line =~ /^gene/i;
    my ($gene, $avc, $wtc, $snf2c, $fc) = split /\t/, $line;
    $gene =~ s/\"//g;
    $gene =~ tr/A-Z/a-z/;
    my $sig = ($fc == 1) ? 0 : 1;
    $h{$gene} = $sig
  }
  return %h; 
  
}

sub ReadFCP
{
  my ($file, $fcsign, $genlist, $colg, $colfc, $colp) = @_;
  
  $colg ||= 1;
  $colfc ||= 2;
  $colp ||= 3;
    
  #print "$file, $fcsign, $genlist, $colg, $colfc, $colp\n";
    
  my %sel = ();
  if(defined $genlist)
  {
    my @selgen = ReadGeneList($genlist);
    %sel = map {$_ => 1} @selgen;
    #print "Selecting ", scalar @selgen, " genes\n";
  }
  #print join(", ", keys %sel), "\n";
  
  my @f = ();
  my @p = ();
  my @genes = ();
  open F, $file or die "Cannot open $file\n";
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    next if $line =~ /^gene/i;
    my @s = split /\t/, $line;
    my $gene = $s[$colg - 1];
    $gene =~ tr/A-Z/a-z/;
    #print "###$gene\n" if $gene =~ /ercc/;
    next if defined $genlist && !defined $sel{$gene};
    
    my $p = $s[$colp - 1];
    my $f = $s[$colfc - 1];
    $p = 1 if $p =~ /NA/;
    $f = 0 if $f =~ /NA/;
    push @p, $p;
    push @f, $f;
    push @genes, $gene;
    #print "$gene $p $f\n";
  }
  die "No genes read.\n" if @f == 0;
  my $lfc = pdl(@f);
  my $p = pdl(@p);
  $lfc *= $fcsign;
  
  return $p, $lfc, \@genes
}


################################################################

=item NormalizeData

  my $nd = NormalizeData($d, $method);
  
  Normalizes data read using ReadExpressionFile: two-dimensional piddle with genes in rows and replicates in columns.
  
  Methods are 'none', 'deseq', 'totalcount', 'tmm', 'fpkm'.

=cut

sub NormalizeData
{
  my ($d, $method, $gff, $genes, $totcnt, $spnorm) = @_;
  
  if(!defined $method || $method eq 'none')
  {
    my ($nrep, $ngene) = dims($d);
    my $f = zeroes($nrep) + 1;
    return $d, $f
  }
  elsif($method eq 'deseq')
  {
    my ($nrep, $ngene) = dims($d);
    my $N = sumover($d>0);
    my $sel = which($N == $nrep);
    my $dd = $d(,$sel);
    my $geom = exp(1/$nrep * sumover(log($dd)))->transpose();
    my $dg = $dd / $geom;
    my ($m, $s, $med) = statsover($dg->transpose());
    my $x = $d / $med;
    #printf "%4.2f ", at($med, $_) for (0 .. nelem($med)-1);print "\n";
    #print "d=$d\nn=$N\ndd=$dd\ngeom=$geom\ndg=$dg\nmed=$med\nx=$x\n\n";
    return $x, $med
  }
  elsif($method eq 'totcountm')
  {
    my ($m) = stats($totcnt);
    $totcnt /= $m;
    my $x = $d / $totcnt;
    return $x, $totcnt
  }
  elsif($method eq 'totcount')
  {
    $totcnt /= 1e6;
    my $x = $d / $totcnt;
    return $x, $totcnt
  }
  elsif($method eq 'loccount')
  {
    my $loccnt = sumover(transpose($d));
    my ($m) = stats($loccnt);
    print "$loccnt\n";
    $loccnt /= $m;
    print "$loccnt\n";
    my $x = $d / $loccnt;
    return $x, $loccnt
  }
  elsif($method eq 'fpkm')
  {
    carp "Need GFF file for FPKM normalization.\n" unless defined $gff && -e $gff;
    my ($nrep, $ngene) = dims($d);
    my %gl = ReadGeneLengths($gff);
    my $length = zeroes($ngene);
    for my $i (0 .. $ngene - 1)
    {
      my $gene = $genes->[$i];
      $gene =~ tr/A-Z/a-z/;
      carp "Missing GFF information for $gene\n" unless defined $gl{$gene};
      $length($i) .= $gl{$gene};
    }
    $length = transpose($length) / 1e3;
    my $totcnt = sumover(transpose($d)) / 1e6;
    my $x1 = $d / $totcnt;
    my $x = $x1 / $length;
    return $x;
  }
  elsif($method eq 'tmm')
  {
    my ($x, $f) = NormalizeRTMM($d)
  }
  elsif($method eq 'tmm_')
  {
    my $Mtrim = 0.3;
    my $Atrim = 0.05;
    
    #$d->badflag(1);$d->badvalue(0);
    my ($nrep, $ngene) = dims($d);
    my $N = sumover($d->transpose());
    my $y = log($d / $N) / log(2);
    my $f = zeroes($nrep) + 1;
    for my $i (1 .. $nrep - 1)
    {
      my $sel = which(($d(0,;-)>0) & ($d($i,;-)>0));
      my $M = ($y(0,$sel) - $y($i,$sel))->flat();
      my $A = (($y(0,$sel) + $y($i,$sel)) / 2)->flat();
      
      my $N1 = at($N, 0);
      my $N2 = at($N, $i);
      my $x1 = $d(0,$sel;-);
      my $x2 = $d($i,$sel;-);
      my $w = ($N1 - $x1)/($N1*$x1) + ($N2 - $x2)/($N2*$x2);
      
      my ($M1, $M2) = trimlimits($M, $Mtrim);
      my ($A1, $A2) = trimlimits($A, $Atrim);
      my $trim = which(($M >= $M1) & ($M <= $M2) & ($A >= $A1) & ($A <= $A2));
      
      my $meanM = sum($w($trim) * $M($trim)) / sum($w($trim));
      my $TMM = 2**$meanM;
      $f($i) .= $TMM;
      
      #print "d=$d\nN=$N\ny=$y\nM=$M\nA=$A\nw=$w\n\n";
      #print "M1=$M1\nM2=$M2\n";
      #print "A1=$A1\nA2=$A2\n";
      #print "meanM=$meanM\ntrim=$trim\n";
      #print "TMM = $TMM\n\n";
      #die;
    }
    $f *= $N;   # this is not in the paper!
    my ($m) = stats($f);
    $f /= $m;           
    my $x = $d / $f;
    #print "$f\n";
    return $x, $f;
  }
  elsif($method eq 'spikein')
  {
    #print "S = $spnorm\n", nelem($spnorm), "\n\n";
    my $x = $d / $spnorm;
    return $x, $spnorm;
  }
  elsif($method eq 'arbitrary')  # do not use in general case!!!
  {
    my $arbcnt = rcols 'arbitrary_norm.txt', 1;
    my ($m) = stats($arbcnt);
    $arbcnt /= $m;
    print "$arbcnt\n";
    my $x = $d / $arbcnt;
    return $x, $arbcnt
  }
  else
    {die "Don't recognize normalization $method.\n"}
}

sub trimlimits
{
  my ($x, $lim) = @_;
  my $s = qsort $x;
  my $n = nelem($s);
  my $nlim = floor($n * $lim + 0.5);
  #print "x=$x\ns=$s\nn=$n\nnlim=$nlim\n";
  my $l1 = at($s, $nlim);
  my $l2 = at($s, $n - $nlim - 1);
    
  #print "l1=$l1\nl2=$l2\n";die;
  return $l1, $l2;
}


sub GenHash
{
  my ($genes, $sel) = @_;
  
  my $n = @$genes;
  $sel = sequence($n) unless defined $sel;
  
  my %h = ();
  for my $i (0 .. nelem($sel) - 1)
  {
    my $idx = at($sel, $i);
    my $gene = $genes->[$idx];
    $gene =~ tr/A-Z/a-z/; 
    $h{$gene} = $i; 
  }
  return %h;
}


sub ReadGeneLengths
{
  my ($gff) = @_;
  
  local *F;
  open F, $gff or die;
  my %h = ();
  while(my $line = <F>)
  {
    chomp $line;
    my ($seq, $src, $feat, $start, $end, $sc, $str, $fr, $attr) = split /\t/, $line;
    next unless $feat eq 'exon';
    my $len = abs($end - $start) + 1;
    my @attr = split(/;/, $attr);
    my $gene = '';
    for my $attr (@attr)
    {
      my ($key, $val) = split " ", $attr;
      $val =~ s/\"//g;
      if($key eq 'gene_id')
        {$gene = $val}
    }
    $gene =~ tr/A-Z/a-z/;
    $h{$gene} = $len;
    #print "$gene: $len\n";
  }
  return %h;
}

sub ReadGFF
{
  my ($gff) = @_;
  
  local *F;
  open F, $gff or die;
  my %gendat = ();
  my (%starts, %ends, %genes);   # arrays for quick gene finding
  while(my $line = <F>)
  {
    chomp $line;
    my ($chr, $src, $feat, $start, $end, $sc, $str, $fr, $attr) = split /\t/, $line;
    my @attr = split(/;/, $attr);
    my $gene_id = '';
    my $gene_name = '';
    for my $attr (@attr)
    {
      my ($key, $val) = split(" ", $attr);
      $val =~ s/\"//g;
      if($key eq 'gene_id')
        {$gene_id = $val}
      if($key eq 'gene_name')
        {$gene_name = $val}
    }
    $gene_id =~ tr/A-Z/a-z/;
    $gene_name =~ tr/A-Z/a-z/;
    push @{$gendat{$gene_id}}, [$feat, $chr, $start, $end, $str, $gene_name];
  }
  
  # find start and end of the entire gene
  for my $gene_id (keys %gendat)
  {
    my $start = 1e16;
    my $end = -1;
    my $chr;
    for my $g (@{$gendat{$gene_id}})
    {
      my ($f, $ch, $st, $en) = @$g;
      $chr = $ch;
      $start = $st if $st < $start;
      $end = $en if $en > $end;
    }
    push @{$starts{$chr}}, $start;
    push @{$ends{$chr}}, $end;
    push @{$genes{$chr}}, $gene_id;
  }
  
  # indexing for quick gene finding
  my %geneidx = ();
  for my $chr (keys %genes)
  {
    my $starts = pdl(@{$starts{$chr}});
    my $ends = pdl(@{$ends{$chr}});
    my $idx = qsorti $starts;
    $geneidx{$chr} = [$starts($idx), $ends($idx), $idx, $genes{$chr}]
  }
   
  return \%gendat, \%geneidx;
}

sub ReadNoGeneFile
{
  my $file = shift;
  
  local *F;
  open F, $file or die "Cannot open file $file\n";
  my @cnt = ();
  while(my $line = <F>)
  {
    chomp $line;
    my @s = split " ", $line;
    shift @s;
    push @cnt, \@s;
  }
  my $cnt = pdl(@cnt);
  my $S = sumover(transpose($cnt));
}

# Simplified: this assumes that there is one-to-one correspondence
# between identifiers and gene names

sub ReadGeneCross
{
  my $file = shift;
  
  local *F;
  open F, $file or die "Cannot open $file\n";
  
  my %g2n = ();
  my %n2g = ();
  while(my $line = <F>)
  {
    chomp $line;
    my ($gene, $name) = split " ", $line;
    $gene =~ tr/A-Z/a-z/;
    $name =~ tr/A-Z/a-z/;
    $g2n{$gene} = $name;
    $n2g{$name} = $gene;
  }
  return \%g2n, \%n2g
}


######################################

sub ReadGOFile
{
  my ($file) = @_;
  
  my %gene2go = ();
  my %go2gene = ();
  my %goname = ();
  
  local *F;
  open F, $file or die "Cannot open $file\n";
  <F>;  # header
  while(my $line = <F>)
  {
    chomp $line;
    my ($gene, $trans, $gname, $go, $goname) = split(/\t/, $line);
    $gene =~ tr/A-Z/a-z/;
    next unless defined $go && $go ne '';
    push @{$gene2go{$gene}}, $go;
    push @{$go2gene{$go}}, $gene;
    $goname{$go} = $goname;
  }
  return \%gene2go, \%go2gene, \%goname
}

################################################################

sub GeneStructure
{
  my ($gene, $gene_data, $quiet) = @_;
  
  die "No gene $gene found in GFF\n" unless defined $gene_data->{$gene};  

  my $start = 1e32;
  my $end = 0;
  my @exons = ();

  my @gg = @{$gene_data->{$gene}};
  my $chromosome;
  my $strand;
  my $gene_name;
  print "$gene\n" unless $quiet;
  for my $g (@gg)
  {
    my ($feat, $chr, $s, $e, $str, $gname) = @$g;
    $chromosome ||= $chr;
    $strand ||= $str;
    $gene_name ||= $gname;
    print "$feat: $chr:$s-$e \($str)\n" unless $quiet;
    $start = $s if ($s < $start);
    $end = $e if ($e > $end);
    push @exons, [$s, $e, $str] if $feat eq 'exon';
  }
  print "$gene_name  $chromosome:$start-$end\($strand)\n" unless $quiet;
  print "\n" unless $quiet;
  return $chromosome, $start, $end, $strand, $gene_name, \@exons
}

################################################################

sub NormalizeRTMM
{
  my $d = shift;
  my ($nrep, $ngene) = dims($d);  
  my $N = sumover($d->transpose());
  
  my $tmp = './normalizeTMM_tmp';
  mkdir "$tmp", 0777 unless -d "$tmp";
  
  my $infile = "$tmp/tmmin$$.dat"; 
  my $outfile = "$tmp/tmmout$$.dat";
  my $Rfile = "$tmp/tmm$$.R";
  unlink $infile, $outfile, $Rfile;

  #my @d = map {$d($_,;-)} (0 .. $nrep-1);  
  #wcols sequence($ngene), @d, $infile;
  local *F;
  open F, ">$infile" or die;
  for my $i (0 .. $ngene - 1)
  {
    my $d = join("\t", list($d(,$i;-)));
    print F "$d\n" 
  }
  close F;

  local *R;
  open R, ">$Rfile" or die "Cannot open temporary file $Rfile\n";
  print R <<EOF;
library("edgeR")
x <- read.delim("$infile", header=F)
group <- factor(rep(1, $nrep))
y <- DGEList(counts=x, group=group)
f <- calcNormFactors(y)
print(f\$samples)
warnings()
EOF
  close R;
  call("R --vanilla < $Rfile >& $outfile", my $stat, 1);
  die "Problem with R call fitdistr\n"if !-e $outfile || $stat;
 
  my @f = ();
  open F, $outfile or die;
  while(my $line = <F>)
    {last if $line =~ /print\(f/}
  <F>;
  while(my $line = <F>)
  {
    chomp $line;
    last if $line =~ /warnings/;
    my ($v, $g, $siz, $f) = split " ", $line;
    push @f, $f;
  }
  close F;
  unlink $infile, $outfile, $Rfile;
  
  my $f = pdl(@f);
  $f *= $N;
  my ($m) = stats($f);
  $f /= $m;           
  my $x = $d / $f;
  return $x, $f
}

################################################################

sub CountFile
{
  my ($experiment, $condition, $replicate, $type) = @_;
  
  my $r = (defined $replicate) ? sprintf("_rep%02d", $replicate) : '';
  my $t = (defined $type) ? "_$type" : '';
  my $file = "$countsdir/${condition}${r}${t}.tsv";
}

################################################################

sub TotCountFile
{
  my ($experiment, $condition, $replicate, $type) = @_;
  
  my $file = "$countsdir/${condition}_readcount.dat";
}

#########################################################

sub SpikeinFile
{
  my ($experiment, $condition, $replicate, $type) = @_;
  
  my $file = "$countsdir/${condition}_spikeinp_normfac.txt";
}

#########################################################
sub OutliersFile
{
  my ($exp, $cond, $rep, $norm, $limit) = @_;
  $limit ||= 5;
  
  my $rr = sprintf "%02d", $rep;
  my $file = "$outlierdir/${cond}_${norm}_rep${rr}_s${limit}.dat";
}

#########################################################


=item GetNormalizedData

  my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel, $f) = GetNormalizedData(%options);

  Read data from expression file, select genes/replicates and normalize.
  
  Returns the following:
  
  $d - two-dimensional data piddle (rows = genes, columns = replicates)
  $m, $s - mean and standard deviation of genes
  $nrep - number of replicates (columns)
  $ngen - number of genes
  $genes - full list (array reference) of genes, before selection
  $name - string to display in figures
  $sel - selection of genes, can be used with $genes array
  $f - normalizing factors
  
  Options:
  
  experiment=>string - not used yet
  condition=>string - condition name (WT, Snf2...)
  type=>string - type of data (raw, fake)
  nonzero=>n - required number of non-zero replicates (-1 = all)
  norm=>string - normalization (deseq, tmm, totalcount, none)
  exclude=>string - exclude a list of replicates from data (comma separated)
  clip=>n - clip n top and n bottom replicates  

=cut

sub GetNormalizedData
{
  my %opt = @_;
  
  $opt{experiment} ||= 'yeast';
  #$opt{condition} ||= 'WT';
  $opt{type} ||= 'raw';
  $opt{nonzero} = 1 unless defined $opt{nonzero};
  $opt{norm} ||= 'deseq';
  $opt{gff} ||= $gfffile;
  
  # special case: levelling normalization is stored in different file
  if($opt{norm} eq 'lev')
  {
    $opt{norm} = 'none';
    $opt{type} = 'lev'
  }

  # include and exclude hashes are {cond}->list, e.g. {WT}->[0,2,5,8,12]

  if($opt{clean})
    {$opt{exclude} = $CleanExclude}
  if($opt{spclean})
    {$opt{exclude} = $SpikeinCleanExclude}
  my %e = HashList($opt{exclude}, 'exclude');
  my %i = HashList($opt{include}, 'include');

  my $name = "$opt{experiment}";
  $name .= " $opt{condition}" if defined $opt{condition}; 
  $name .= " Rep$opt{replicate}" if defined $opt{replicate};
  $name .= ", ex $opt{exclude}" if $opt{exclude};
  $name .= ", in $opt{include}" if $opt{include};

  # nasty way of incorporating same-condition bootstrap
  # split condition into two random subsets, and pass the selections
  # through include lists %i
  my @myconds = @conds;
  if(defined $opt{randrep2})
  {
    die "Need specified condition for randrep2\n" unless defined $opt{condition};
    @myconds = ("$opt{condition}.1", "$opt{condition}.2");
    my $file = CountFile($opt{experiment}, $opt{condition}, $opt{replicate}, $opt{type});
    my $nrep = Nrep($file);
    my ($xsel1, $xsel2) = SplitCondition($e{$opt{condition}}{exclude}, $nrep);
    #print "S1=", qsort($xsel1), "\nS2=", qsort($xsel2), "\n";
    my $nsplit = nelem($xsel1);
    
    my ($c1, $c2) = @myconds;
    $i{$c1}{include} = [list($xsel1+1)];
    $i{$c2}{include} = [list($xsel2+1)];
    $opt{randrep} = $opt{randrep2};
    die "-randrep2=$opt{randrep2} is too large. Can only select $nsplit reps.\n" if $opt{randrep2} > $nsplit;
  }

  # read both conditions, so we can get normalization
  my %dat = ();
  my %repsel = ();
  my @genes = ();
  my %totcnt = ();
  my %spnorm = ();
  my $div = -1;
  for my $cond (@myconds)
  {
    # $cond - indexing condition, $fcond - true condition (for same-cond bootstrap only)
    my $tcond = (defined $opt{randrep2}) ? $opt{condition} : $cond;
    my $file = CountFile($opt{experiment}, $tcond, $opt{replicate}, $opt{type});
    my $cntfile = TotCountFile($opt{experiment}, $tcond);
    my $spfile = SpikeinFile($opt{experiment}, $tcond);
    print "Reading data from $file..." unless $opt{quiet};
    my ($d, $m, $s, $sel, $repsel) = ReadExpressionFile($file, nonzero=>0, %{$e{$cond}}, %{$i{$cond}}, randrep=>$opt{randrep});
    my @g = ReadGeneList($file);
    my $cnt = ReadTotCountFile($cntfile, repsel=>$repsel);
    my $spn = ReadSpikeinFile($spfile, repsel=>$repsel);
    print "done\n" unless $opt{quiet};
    $dat{$cond} = $d;
    $repsel{$cond} = $repsel;
    @genes = @g;        # should be the same, not checking!
    $totcnt{$cond} = $cnt;
    $spnorm{$cond} = $spn;
    my ($nrep, $ngen) = dims($d);
    $div = $nrep if $div == -1;   # where to divide conditions
  }

  # merge data for normalization
  my ($c1, $c2) = @myconds;
  my $alldat = append($dat{$c1}, $dat{$c2});
  my $allcnt = append($totcnt{$c1}, $totcnt{$c2});
  my $allspnorm = append($spnorm{$c1}, $spnorm{$c2});
  my ($nrepall, $ngenall) = dims($alldat);
  print "$nrepall x $ngenall in total\n" unless $opt{quiet};
  
  # normalize all data  
  print "Normalizing data by $opt{norm}..." unless $opt{norm} eq 'none' || $opt{quiet};
  my ($nalldat, $allf) = NormalizeData($alldat, $opt{norm}, $opt{gff}, \@genes, $allcnt, $allspnorm);
  print " done\n" unless $opt{norm} eq 'none' || $opt{quiet};
  
  # separate conditions
  my %D = ();
  
  for my $cond (@myconds)
  {
    my ($d, $f);
    if($cond eq $c1)
    {
      $d = $nalldat(0:$div-1,);
      $f = $allf(0:$div-1) unless $opt{norm} eq 'fpkm';
    }
    elsif($cond eq $c2)
    {
      $d = $nalldat($div:-1,);
      $f = $allf($div:-1) unless $opt{norm} eq 'fpkm';
    }
    my $repsel = $repsel{$cond};
    my ($nrep, $ngen) = dims($d);
    print "### $cond: $nrep x $ngen\n";
    print "### Selection: ", qsort($repsel)+1, "\n";
    
    # selection on number of zeroes
    my $sel = sequence($ngen);
    if(defined $opt{nonzero})
    {
      my $nz = ($opt{nonzero} < 0) ? $nrep : $opt{nonzero};
      my $N = sumover($d>0);
      $sel = which($N >= $nz);
      $d = $d(,$sel);
      $ngen = nelem($sel);
    }
    
    if($opt{clip})
    {
      my ($nrep, $ngen) = dims($d);
      my $t1 = $opt{clip};
      my $t2 = $nrep - $opt{clip} - 1;
      if($t2 > $t1)
      {
        $d = qsort $d;
        $d = $d($t1:$t2);
        $name .= ", clip $opt{clip}";
      }
    }
    
    my ($m, $s) = statsover($d);
    ($nrep, $ngen) = dims($d);
    print "Selected $ngen genes in $nrep replicates in $cond.\n" unless $opt{quiet};
    $D{$cond} = [$d, $m, $s, $nrep, $ngen, \@genes, $name, $sel, $f, $repsel];
  }
  
  if(defined $opt{condition} && !defined $opt{randrep2})
    {return @{$D{$opt{condition}}}}
  else
    {return %D} 
}


#
# $list can be "1,3,5:" or ":1,2,3" or "1:3,4" 
#
sub HashList
{
  my ($list, $opt) = @_;
  
  my %h = ();
  if(defined $list)
  {
    die "Replicate list $list must contain colon\n" unless $list =~ /:/;
    my @s = split(/:/, $list);
    for my $i (0 .. 1)
    {
      next unless defined $s[$i];
      my @a = split(/,/, $s[$i]);
      next if @a == 0;
      $h{$conds[$i]}{$opt} = \@a
    }
  }
  return %h
}


sub HashList_
{
  my ($list, $opt, $condition) = @_;
  
  my %h = ();
  if(defined $list)
  {
    my @s = split(/:/, $list);
    if(@s == 1)   # only current condition
    {
      my @a = split(/,/, $list);
      $h{$condition}{$opt} = \@a
    }
    elsif(@s == 2)  # both conditions specified
    {
      for my $i (0 .. 1)
      {
        my @a = split(/,/, $s[$i]);
        $h{$conds[$i]}{$opt} = \@a
      }
    }
  }
  return %h
}

sub Nrep
{
  my $file = shift;
  
  local *F;
  open F, $file or die "Cannot open file $file\n";
  my $line = <F>;
  chomp $line;
  my @s = split /\t/, $line;
  my $nrep = scalar(@s) - 1;
  return $nrep
}

sub SplitCondition
{
  my ($exclude, $nrep) = @_;
  
  my $xsel;
  if(defined $exclude)
  {
    my $ex = pdl($exclude) - 1;
    my $all = sequence($nrep);
    $xsel = setops($all, 'XOR', $ex);
  }
  else
    {$xsel = sequence($nrep)}
  my $n = nelem($xsel);
  my $n2 = floor($n / 2);

  my $r = qsorti random $n;
  my $r1 = $r(0:$n2-1);  # random selection
  my $r2 = $r($n2:$n-1);
  my $xsel1 = $xsel($r1);
  my $xsel2 = $xsel($r2);
  return $xsel1, $xsel2;
}

################################


=item PlotDistWithFits

  PlotDistWithFits($win, $pan, $d, $lab, %opt);

  Fit data with theoretical distributions and plot histogram/KDE of data and best-fitting distributions.
  
  $win - window object (see PLGraphs)
  $pan - panel number (1, 2, 3, ...)
  $d - data (one-dimensional piddle)
  $lab - label in the figure
  %opt - options:
    charsize - character size
    hist - plot histogram
    kde - plot KDE
    onebin - histogram in bins of 1
    bestbin - select best binning automatically
    max - maximum in y axis 

=cut


sub PlotDistWithFits
{
  my ($win, $pan, $d, $lab, %opt) = @_;
  
  $lab ||= '';
  $opt{charsize} ||= 0.6;
  my $cs = $opt{charsize};
      
  my $mx = max($d)+1;
  my $p1 = -1.5;
  my $p2 = $mx + 1.5;
  my $nx = $mx+3;
  
  $opt{onebin} ||= $mx < 100;
  undef $opt{onebin} if $opt{bestbin};
  
  my ($x, $h);
  if($opt{onebin})
    {($x, $h) = BuildHistogram($d, $nx, $p1, $p2, 1)}
  else
    {($x, $h) = BuildHistogram($d, undef, undef, undef, 1, 1)}
  my ($xk, $kde) = KDE($d, undef, 10);
  
  my ($M, $S) = stats($d);
  my $nozeroes = (nelem($d->where($d==0)) == 0);
  $opt{kde} ||= ($M > 10);
  undef $opt{kde} if $opt{hist};
  
  my $sk = Skewness($d);
  my $cv = $S / $M;
  
  print "Fitting...";
  my ($nb_loglik, %nb_pars) = FitDistr(floor($d+0.5), 'negative binomial');
  my ($ln_loglik, %ln_pars) = FitDistr($d->where($d>0), 'lognormal') if $nozeroes;
  my ($nm_loglik, %nm_pars) = FitDistr($d, 'normal');
  #my ($ga_loglik, %ga_pars) = FitDistr($d->where($d>0), 'gamma');
  print "done\n";

  my $loglikrat = 2 * (-$nb_loglik + $ln_loglik)  if $nozeroes;

  print "Calculating distributions...";
  my ($nb_M, $nb_S) = GetRMS('negative binomial', %nb_pars);
  my ($ln_M, $ln_S) = GetRMS('lognormal', %ln_pars) if $nozeroes;
  my ($nm_M, $nm_S) = GetRMS('normal', %nm_pars);
  #my ($ga_M, $ga_S) = GetRMS('gamma', %ga_pars);
  print "done \n";
  
  print "NB: M=$nb_M S=$nb_S\n";
  print "LN: M=$ln_M S=$ln_S\n"  if $nozeroes;
  print "NM: M=$nm_M S=$nm_S\n";
  
  my ($nb_L, $nb_sic) = lnL(floor($d+0.5), 'negative binomial', %nb_pars);
  my ($ln_L, $ln_sic) = lnL($d->where($d>0), 'lognormal', %ln_pars) if $nozeroes;
  my ($nm_L, $nm_sic) = lnL($d, 'normal', %nm_pars);
  #my ($ga_L, $ga_sic) = lnL($d->where($d>0), 'gamma', %ga_pars);

  print "NB: L1 = $nb_loglik  L2 = $nb_L\n";
  print "LN: L1 = $ln_loglik  L2 = $ln_L\n" if $nozeroes;
  print "NM: L1 = $nm_loglik  L2 = $nm_L\n";
  
  #my ($mn, $sn) = mcfit($d, \&nbinom, 1, 100);
  #my ($ml, $sl) = mcfit($d, \&lnorm, 0, 500);

  print "Testing...";
  my ($Pn, $Pl, $Pm) = ('N/A', 'N/A', 'N/A');
  my $Xf;
  if($nozeroes)
    {$Pl = NormalityTest(log10($d))}
  $Pn = NBBootstrapTest('m', $d, 5, 1e3) unless $opt{simple};
  $Pm = NormalityTest($d);
    
  $Pl = sprintf "%.3g", $Pl unless $Pl eq 'N/A';
  $Pn = sprintf "%.3g", $Pn unless $Pn eq 'N/A';
  $Pm = sprintf "%.3g", $Pm unless $Pm eq 'N/A';
  $Pn = $opt{Pn} if defined $opt{Pn};
 
  #my ($Dn, $Pn) = KSone($d, \&nbinomms, $nb_M, $nb_S); # cannot do discrete!
  #my ($Dl, $Pl) = KSone($d, \&lnorm, $ln_M, $ln_S) if $nozeroes;
  #my ($chin, $Pn) = ChiSquareDiscrete($d, \&fnbinomms, $nb_M, $nb_S);
  print "done\n";
  
  #print "NB: D = $Dn, P = $Pn\n";
  #print "LN: D = $Dl, P = $Pl\n";
  
  print "Calculating models...";
  my ($x_nb, $y_nb) = CreateDiscreteDist(\&fnbinomms, $nb_M, $nb_S) if $nb_M;
  my ($x_ln, $y_ln) = CreateDist(\&flnorm, $ln_M, $ln_S) if $nozeroes && $ln_M;
  my ($x_nm, $y_nm) = CreateDist(\&fgauss, $nm_M, $nm_S) if $nm_M;
  print " done\n";
  
  my $Q = $y_nb;
  $Q = append($Q, $y_ln) if $nozeroes;
  $Q = append($Q, $y_nm);
  $Q = append($Q, ($opt{kde}) ? $kde : $h);
  print nelem($Q), "\n";
  my $max = max($Q) * 1.2;
  $max = 1.1*$max unless $opt{simple};
  $max = $opt{max} if defined $opt{max};
  my $min = 0;
  my $minx = -0.5;
  my $maxx = max($x)*1.2;
  
  my ($xf, $xmul) = pfac($maxx);
  my ($yf, $ymul) = pfac($max, 2);
  
  $minx /= $xf;
  $maxx /= $xf;
  $max /= $yf;
  
  my ($xcloud, $ycloud);
  if($opt{discrete})
  {
    my @cnt = (0) x max($d);
    for my $i (0 .. nelem($d) - 1)
    {
      my $dd = at($d, $i);
      $cnt[$dd]++
    }
    #print "CNT=", pdl(@cnt), "\n";
    
    my @xc = (); my @yc = ();
    for my $i (0 .. @cnt-1)
    {
      my $c = $cnt[$i];
      if($c > 0)
      {
        for my $k (0 .. $c-1)
        {
          push @xc, $i;
          push @yc, $k - ($c-1)/2;
        }  
      }
    }
    $xcloud = pdl(@xc) / $xf;
    $ycloud = pdl(@yc);
    my $cmax = max($ycloud);
    $ycloud /= $cmax if  $cmax > 0;
    #print "$d\n";
    #wcols $xcloud, $ycloud;
    #die;
  }
  else
  {
    $ycloud = grandom($d) * 0.3;
    $ycloud->where($ycloud > 1) .= 0;
    $ycloud->where($ycloud < -1) .= 0;
    $xcloud = $d / $xf;
  }


  my $xpan = ($pan-1) % $PL_nx + 1;
  my $ypan = int(($pan-1) / $PL_nx) + 1;
  my ($x1, $x2, $y1, $y2) = PanelCoord($xpan, $ypan);
  my $ysplit = $y1 + 0.1 * ($y2 - $y1);

  my $xpow = ($xmul == 0) ? '' : " (x1e$xmul)";
  my ($xtick, $nxsub) = (max($x) > 20) ? (0, 0) : (2,2);
   
  $win->xyplot(pdl(-1),pdl(-2),
    VIEWPORT=>[$x1, $x2, $y1, $ysplit],
    BOX => [$minx, $maxx, -1.1, 1.1],
    XLAB => "Count$xpow", YLAB => '',
    XBOX => 'BINST', YBOX =>'B',
    MINTICKSIZE=>0.5, MAJTICKSIZE=>0.5,
    XTICK=>$xtick, NXSUB=>$nxsub,
    CHARSIZE=>$cs
   );

  #LinePlot($win, pdl(-1, $maxx), pdl(0, 0));
  PointPlot($win, $xcloud, $ycloud, SYMBOLSIZE=>0.4);
  
  if($opt{outliers})
  {
    my ($clean_dat, $iup, $idown, $M, $S) = RejectOutliersSigma($d/$xf);
    vline($win, $M, 1, 'RED');
    vline($win, $M+5*$S, 1, 'GREY');
    vline($win, $M-5*$S, 1, 'GREY');
  }
  
  sub vline {my ($w,$x, $s, $c) = @_; LinePlot($w, pdl($x,$x), pdl(-2,2), COLOR=>$c, LINESTYLE=>$s); print "X=$x\n"}
  
  $win->xyplot(pdl(-1),pdl(-2),
    VIEWPORT=>[$x1, $x2, $ysplit, $y2],
    BOX => [$minx, $maxx, $min, $max],
    %PL_empty
   );


  if($opt{kde})
    {OverPlot($win, $xk/$xf, $kde/$yf, PLOTTYPE=>'LINE', COLOR=>[110,200,110], LINEWIDTH=>3)}
  else
  {
    BarPlot($win, $x/$xf, $h/$yf, colour=>[215,255,215], y0=>0);
    HistPlot($win, $x/$xf, $h/$yf, colour=>[110,200,110], linewidth=>3);
  }
  #PlotDist(\&fgauss, $M, $S, 'RED');
  #PlotDist(\&fpoisson, $M, $S, 'GREEN');
  #PlotDDist($win, \&fnbinomms, $nb_M, $nb_S, 'BLUE');
  #PlotDist($win, \&flnorm, $ln_M, $ln_S, 'RED') if $nozeroes;
  #PlotDist(\&fgammams, $ga_M, $ga_S, 'GOLD2') if $nozeroes;
  
  HistPlot($win, $x_nb/$xf, $y_nb/$yf, colour=>'BLUE', linewidth=>3) if $nb_M;
  LinePlot($win, $x_ln/$xf, $y_ln/$yf, COLOR=>'RED', LINEWIDTH=>3) if $nozeroes && $ln_M;
  LinePlot($win, $x_nm/$xf, $y_nm/$yf, COLOR=>'GREY', LINEWIDTH=>3) if  $nm_M;
  
  my $ypow = ($ymul == 0) ? '' : " (x1e$ymul)";
  $win->xyplot(pdl(-1),pdl(-2),
    VIEWPORT=>[$x1, $x2, $ysplit, $y2],
    BOX => [$minx, $maxx, $min, $max],
    YLAB => "Frequency$ypow",
    CHARSIZE => $cs,
    XBOX => 'B', YBOX=>'BNISTF',
    MINTICKSIZE=>0.5, MAJTICKSIZE=>0.5
   );

  
  my $yc = 0;
  my $dy = 2.1;
  TopText($win, $lab, CHARSIZE=>1*$cs, y=>$yc, COLOR=>'BLACK');$yc-=$dy;
  #TopText($win, sprintf("M=%.3g S=%.3g", $M, $S), CHARSIZE=>$cs, y=>$yc); $yc-=$dy;
  #TopText($win, sprintf("CV=%.3g G=%.3g N=%d", $cv, $sk, $N), CHARSIZE=>$cs, y=>0);
  #TopText($win, sprintf("Dn = %.3g  Dl = %.3g", $Dn, $Dl), y=>-3, CHARSIZE=>0.4);
  if($opt{withprobs}) {TopText($win, "PNM = $Pm PNB = $Pn PLN = $Pl ", CHARSIZE=>0.7*$cs, y=>$yc); $yc-=$dy;}
  

  #TopText($win, sprintf("SIC#dNB#u = %.3g  SIC#dLN#u = %.3g  D=%.3g", $nb_sic, $ln_sic, $loglikrat), y=>$yc, CHARSIZE=>0.8*$cs);
  unless($opt{simple})
  {
    TopText($win, sprintf("NB: P = %s  SIC = %.3g", $Pn, $nb_sic), CHARSIZE=>$cs, y=>$yc); $yc-=$dy;
    TopText($win, sprintf("LN: P = %s  SIC = %.3g", $Pl, $ln_sic), CHARSIZE=>$cs, y=>$yc); $yc-=$dy;
  }
}

################################

sub pfac
{
  my ($x, $lim) = @_;
  $lim ||= 4;
  
  my $l = log($x)/log(10);
  
  my $fac = 0;
  if($l < -$lim || $l > $lim)
  {
    $fac = floor($l)
  }
  return 10**$fac, $fac 
}

################################

sub PlotDist
{
  my ($win, $f, $M, $S, $col) = @_;
  
  my ($x, $y) = CreateDist($f, $M, $S);
  OverPlot($win, $x, $y*1000, PLOTTYPE=>'LINE', COLOR=>$col, LINEWIDTH=>3); 
}

################################

sub PlotDDist
{
  my ($win, $f, $M, $S, $col) = @_;
  
  my ($x, $y) = CreateDiscreteDist($f, $M, $S);
  HistPlot($win, $x, $y*1000, colour=>$col, linewidth=>3); 
}

#########################################################

sub CreateDist
{
  my ($f, $M, $S, $n, $ns, $shift) = @_;

  $n ||= 1000;
  $ns ||= 10;
  $shift = 1 unless defined $shift;
  
  my $x = ($M + $ns*$S) * (sequence($n) + $shift)/$n;
  $x = $x->where($x>0);
  my $y = zeroes($x);
  for my $i (0 .. $n - 1)
    {$y($i) .= &$f(at($x, $i), $M, $S)}
  
  return $x, $y
}

#########################################################

sub CreateDiscreteDist
{
  my ($f, $M, $S, $cumul) = @_;

  #$n ||= int($M + 10*$S);
  #my $gap = 1;
  #if($n > 1000)
  #  {$gap = int($n / 1000)}
   
  #print "n=$n\ngap=$gap\n";
  #my $x = sequence($n / $gap) * $gap;
  #print "$x\n";
  
  my $x = CreateDiscreteGrid($f, $M, $S, $cumul);
  #wcols "$x\n";
  
  my $y = zeroes($x);
  for my $i (0 .. nelem($x) - 1)
    {$y($i) .= &$f(at($x, $i), $M, $S)}
  
  #wcols $x, $y;
  
  return $x, $y
}

#########################################################

sub RejectOutliersSigma
{
  my ($x, $Z, $trim) = @_;
  $trim = 3 unless defined $trim;
  
  my $q = qsort $x;
  my $n = nelem($x);
  
  carp "trim too large\n" if $trim > $n - $trim - 1;
  my ($m, $s) = stats($q($trim:$n-$trim-1));
  my $up = which($x - $m >= $s * $Z); 
  my $down = which($m - $x >= $s * $Z);
  my $in = setops(sequence($n), 'XOR', append($up, $down));
  my $xx = $x($in)->copy();
  
  return ($xx, $up, $down, $m, $s);
}

#########################################################

sub RejectOutliersChauvenet
{
  my ($x) = @_;
  
  my $iq = qsorti $x;
  my $q = $x($iq);
  my $n = nelem($x);
  my ($m, $s) = stats($q);
 
  my @down = ();
  my $p;
  my $i = 0;
  while($i <= $n - 1)
  {
    my $z = abs(at($q, $i) - $m) / $s;
    $p = 2 * (1 - pnorm($z));
    #print "dn ", at($q, $i), "   $z   $p\n";
    if($n * $p < 0.5)
      {push @down, at($iq, $i)}
    else
      {last}
    $i++
  }
  
  my @up = ();
  $i = $n - 1;
  while($i >= 0)
  {
    my $z = abs(at($q, $i) - $m) / $s;
    $p = 2 * (1 - pnorm($z));
    #print "up ", at($q, $i), "   $z   $p\n";
    if($n * $p < 0.5)
      {push @up, at($iq, $i)}
    else
      {last}
    $i--
  }
  
  my $test = sequence(1);
  my $empty = $test->where($test > 100);
  
  my $down = (@down > 0) ? pdl(@down) : $empty;
  my $up = (@up > 0) ? pdl(@up) : $empty;
  my $in = setops(sequence($n), 'XOR', append($up, $down));
  my $xx = $x($in)->copy();
  
  return ($xx, $up, $down);
}

#########################################################

sub NormalityFromData
{
  my ($p, $log, $outliers, $sigma) = @_;
  
  $sigma ||= 5;
  
  my ($M, $N) = dims($p);
  my @P = ();
  for my $i (0 .. $N - 1)
  {
    Percent1($i/$N) if $i % 100 == 0;
    my $d = $p(,$i;-);
    #my $n = nelem($d);
    $d = log10($d) if $log;
    
    if(defined $outliers)
    {
      if($outliers eq 'chauvenet')
        {($d) = RejectOutliersChauvenet($d)}
      elsif($outliers eq 'sigma')
        {($d) = RejectOutliersSigma($d, $sigma)}
    }
    
    my $P = NormalityTest($d);
      
    push @P, $P;
  }
  print "\n";
  my $P = pdl(@P);

  return $P;  
}

#########################################################

sub PoissonFromData
{
  my ($dat, $type, $quiet) = @_;
  
  my ($M, $N) = dims($dat);
  my @P = ();
  for my $i (0 .. $N - 1)
  {
    unless($quiet)
      {Percent1($i/$N) if $i % 100 == 0}
    my $d = $dat(,$i;-);
    my $n = nelem($d);
    
    my $p = PoissonTest($d, $type);
      
    push @P, $p;
  }
  print "\b\b\b\b\b\b\b done    \n" unless $quiet;
  return pdl(@P);  
}

#########################################################

sub NBHGFromData
{
  my ($dat, $devel) = @_;
  
  my ($M, $N) = dims($dat);
  my @P = ();
  for my $i (0 .. $N - 1)
  {
    Percent1($i/$N) if $i % 100 == 0;
    my $d = $dat(,$i;-);
    my $n = nelem($d);
    
    my ($Xf, $p) = HinzGurland($d, $devel);
    $p = 0 if isnan($p);
    
    push @P, $p;
  }
  print "\b\b\b\b\b\b\b done    \n";
  return pdl(@P);  
}

sub isnan { ! defined( $_[0] <=> 9**9**9 ) }

#################################################

sub NBKSFromData
{
  my ($dat) = @_;
  
  my ($M, $N) = dims($dat);
  my @P = ();
  for my $i (0 .. $N - 1)
  {
    Percent2($i/$N);# if $i % 100 == 0;
    my $d = $dat(,$i;-);
    my $n = nelem($d);
    
    my $R = '/sw/opt/R/2.15.1/bin/R';
    my ($D, $p) = NBKSTest($d, undef, $R);
    push @P, $p;
  }
  print "\b\b\b\b\b\b\b done    \n";
  return pdl(@P);  
}

#################################################

=item B<ReadTestfileInfo>

  my %info = ReadTestfileInfo($file)
  
Reads DE tool information from a metafile $file.

=cut

sub ReadTestfileInfo
{
  my ($file) = @_;
  
  
  local *F;
  open F, $file or die "Cannot open file $file\n";
  
  my %h = ();
  my %t = ();
  my $test;
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    next if $line =~ /^\s*$/;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    if($line =~ /TEST\s+(\S+)/)
    {
      if(defined $test)
      {
        defaults(\%t, $test, $file);
        $h{$test} = {%t};
        #print "$test ### ", join(", ", map {"$_ = $t{$_}"} keys %t), "\n";
        %t = ();
        undef $test;
      }
      $test = $1
    }
    elsif($line =~ /(\S+)\s*=\s*(.+)/)
    {
      $t{$1} = $2;
    }
  }
  defaults(\%t, $test, $file);
  $h{$test} = \%t;
  
  return %h
}

sub defaults
{
  my ($t, $test, $file) = @_;
  my @compulsory = qw(name fullfile);
  my %default = (
    repfile => '-',
    sp_repfile => '-',
    adjustment => 'none',
    'FC-sig' => 1
  );
  for my $c (@compulsory)
    {die "'$c' not defined for test $test in $file\n" unless defined $t->{$c}}
  for my $k (keys %default)
  {
    my $v = $default{$k};
    $t->{$k} = $v unless defined $t->{$k}
  }
  $_[0] = $t
}


sub ReadTestfileInfoOld
{
  my ($file) = @_;
  
  local *F;
  open F, $file or die "Cannot open file $file\n";
  
  my %h = ();
  while(my $line = <F>)
  {
    next if $line =~ /^#/;
    chomp $line;
    my ($id, $name, $bsform, $sp_bsform, $fullfile, $cor, $fcsig) = split /\t/, $line;
    $h{$id} = [$name, $bsform, $sp_bsform, $fullfile, $cor, $fcsig];
  }
  return %h
}

#################################################

sub MulticorLimit
{
  my ($p, $current, $new, $alpha) = @_;
  
  $alpha ||= 0.05;
  
  my $c;
  if($current eq 'none')
    {$c = $new}
  elsif($current eq $new)
    {$c = 'none'}   # do not correct already corrected
  elsif($current =~ /lim(.+)/)  # explicit limit 
  {
    $c = 'none';
    $alpha = $1;
  }
  else
    {die "Cannot correct with $new as it is already corrected with $current\n"}
    
  my ($lim, $ns);
  if($c eq 'hb')
    {($lim, $ns) = HolmBonferroniLimit($p, $alpha)}
  elsif($c eq 'bh')
    {($lim, $ns) = BenjaminiHochbergLimit($p, $alpha)}
  elsif($c eq 'none')
  {
    $lim = $alpha;
    $ns = nelem($p->where($p<=$alpha));
  }
  else
    {die "Unrecognized correction $c\n"}
  #print "Corrected $current to $new by $c at $alpha\n";
  
  return $lim unless wantarray;
  return $lim, $ns
}

####################################

sub rank2
{
  my ($p, $fc) = @_;
  my $fc1 = 1 / (1 + abs($fc));
  my $d = transpose(pdl($p, $fc1));
  my $idx = qsortveci($d);
  my $s = sequence($idx) + 1;
  my $rank = zeroes($s);
  $rank($idx) .= $s;
  return $rank
}

####################################


sub PowerTest
{
  my ($p1, $p2, $fc1, $fc2, $multicor, $curcor1, $curcor2, $fclim, $alpha, $oneside) = @_;
  
  $alpha ||= 0.05;
  
  # one-sided test is for 0 <= fc1 <fclim
  # two-sided test is the other way around: for fc1 >=  fclim
  # these tests are used for different things
  my $sel;
  if($oneside)
  {
    if($fclim > 0)
      {$sel = which(($fc1 < $fclim) & ($fc1 >= 0))}
    else
      {$sel = which(($fc1 > $fclim) & ($fc1 <= 0))}
  }
  else
    {$sel = which(abs($fc1) >= $fclim)}
    
  return 0,0,0,0 if nelem($sel) == 0;
   
  $p1 = $p1($sel);
  $p2 = $p2($sel);

  my $l1 = MulticorLimit($p1, $curcor1, $multicor, $alpha);
  my $l2 = MulticorLimit($p2, $curcor2, $multicor, $alpha);
  
  my $s1 = ($curcor1 =~ /llim/) ? ($p1 >= $l1) : ($p1 <= $l1);
  my $s2 = ($curcor2 =~ /llim/) ? ($p2 >= $l2) : ($p2 <= $l2);
  
  my $TP = nelem(which(($s1 == 1) & ($s2 == 1)));
  my $TN = nelem(which(($s1 == 0) & ($s2 == 0)));
  my $FP = nelem(which(($s1 == 0) & ($s2 == 1)));
  my $FN = nelem(which(($s1 == 1) & ($s2 == 0)));
  
  return $TP, $TN, $FP, $FN, $s1, $s2
}

sub SpikeinPowerTest
{
  my ($true, $p, $multicor, $curcor, $alpha) = @_;
  
  $alpha ||= 0.05;
  
  my $l = MulticorLimit($p, $curcor, $multicor, $alpha);
  my $s = ($p <= $l);
  #wcols $true, $s, $p;
  
  my $TP = nelem(which(($s == 1) & ($true == 1)));
  my $TN = nelem(which(($s == 0) & ($true == 0)));
  my $FP = nelem(which(($s == 1) & ($true == 0)));
  my $FN = nelem(which(($s == 0) & ($true == 1)));
  
  #print "TP=$TP TN=$TN FP=$FP FN=$FN\n";
  
  return $TP, $TN, $FP, $FN
}
  
1;