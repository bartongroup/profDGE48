#!/sw/bin/perl

=head1 NAME

B<plot_gene_examples.pl>

=head1 DESCRIPTION

Plots a figure with expression distribution for several genes. Best-fitting log-normal and negative binomial distributions are shown.

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

#require 'mylib.pl';

$| = 1;
Stamp();

#my @conds = qw(WT Snf2);

$genes = 'yar064w,ybr019c,ybr011c,ycr050c,yer103w,ylr129w';
my ($kde, $hist, $onebin, $bestbin, $max);
my $simple = 1;
GetMyOptions(
  #'max=f' => \$max,
  #kde => \$kde,
  #hist => \$hist,
  #onebin => \$onebin,
  #bestbin => \$bestbin,
  #simple => \$simple
);
die "Need -genes\n" unless defined $genes;

ReadDefs();

###############################################

my %g2i = ();
my %data = GetNormalizedData(
    replicate=>$rep,
    exclude=>$exclude,
    nonzero=>$nonzero,
    norm=>$norm,
    type=>$type,
    clean=>$clean,
    spclean=>$spclean
);
for my $c (@conds)
{
  my %g = GenHash($data{$c}[5], $data{$c}[7]);
  $g2i{$c} = \%g;
  # $d, $m, $s, $nrep, $ngen, $genes, $name, $sel 
}

######################################################

my @genes = split(/,/, $genes);
my $N = @genes;

my %nb = ReadNBTest($cond);

$PL_xmin = 0.1;
$PL_ymin = 0.3;
$PL_nx = 3;
$PL_ny = 2;
$PL_deltax = 0.06;
$PL_deltay = 0.06;
my $cs = ($simple) ? 0.5 : 0.45;
my $w = NewWindow($psfile);

#my $pan = $N**2 - $N + 1;
my $pan = 1;
for my $gene (@genes)
{
  $gene =~ tr/A-Z/a-z/;
  die "Missing data for $gene\n" unless defined $g2i{$cond}{$gene};
  
  my $GENE = $gene;
  $GENE =~ tr/a-z/A-Z/;
  my $title = "(" . chr(96+$pan) . ") " . $GENE;
  
  my ($d, $m, $s, $nrep, $ngen, $gns, $sname, $sel) = @{$data{$cond}};
  if(defined $g2i{$cond}{$gene})
  {
    my $i = $g2i{$cond}{$gene};
    my $dat = $d(,$i)->flat();
    my $Pn = sprintf "%.2g", $nb{$gene};
    PlotDistWithFits($w, $pan, $dat, $title, 
      kde=>$kde,
      hist=>$hist,
      onebin=>$onebin,
      bestbin=>$bestbin,
      charsize=>$cs,
      simple=>$simple,
      withprobs=>1,
      Pn=>$Pn,
      discrete=>1
    );
  }
  $pan++
}

$w->close();


sub ReadNBTest
{
  my ($cond) = @_;
 
  my $file = "./${cond}_nbtest_m_clean.dat";
  local *F;
  open F, $file or die;
  my %h = ();
  while(my $line = <F>)
  {
    chomp $line;
    my ($gene, $m, $s, $p) = split /\t/, $line;
    $gene =~ tr/A-Z/a-z/;
    $h{$gene} = $p;
  }
  return %h
}


=head1 SYNOPSIS

  plot_gene_examples.pl -clean -psfile=figure.ps

=head1 OPTIONS

=over 4

=item B<-clean>

Use clean data replicates, as specified in C<defs.dat> file.

=item B<-psfile>=I<name>

Name of the postsript file to redirect output to.

=back
