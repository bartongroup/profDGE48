#!/sw/bin/perl

=head1 NAME

B<one_nb_test.pl>

=head1 DESCRIPTION

This scrip is a part of C<grid_launcher_nbtest.pl>. It is not supposed to be run directly.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;


use PDL;
use PDL::NiceSlice;

use GRNASeq;
use GetOpt;

use CompBio::Tools;
use CompBio::Distribution;

$| = 1;

my $created = "Created with $0 " . join(" ", @ARGV) . "\n";

ReadDefs();

my ($ncrit, $maxn);
my $batchsize = 30;
my $stat = 'm';
my $batch = 0;
my $outfile;
my ($help, $man);
GetMyOptions(
  'batch=i' => \$batch,
  'stat=s' => \$stat,
  'outfile=s' => \$outfile,
  'batchsize=i' => \$batchsize,
  'ncrit=i' => \$ncrit,
  'maxn=i' => \$maxn
);

print "$created\n";
print "Batch $batch\n";
print "Outfile $outfile\n";
print "ncrit $ncrit\n" if defined $ncrit;
print "maxn $maxn\n" if defined $maxn;

my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel) = GetNormalizedData(
  condition=>$cond,
  replicate=>$rep,
  exclude=>$exclude,
  nonzero=>1,
  norm=>$norm,
  type=>$type,
  clean=>$clean,
  spclean=>$spclean
);


####################################

#$batch++ if $batch == 0;

my $n1 = ($batch - 1) * $batchsize;
my $n2 = $n1 + $batchsize - 1;
$n2 = $ngen - 1 if $n2 > $ngen - 1;
exit if $n2 < $n1;

print "Processing batch $batch, $n1 .. $n2 \(n = $ngen)\n";

local *F;
open F, ">$outfile" or die;


for my $i ($n1 .. $n2)
{
  my $gene = $genes->[at($sel, $i)];
  my @out = ($gene);
  my $x = $d(,$i;-);
  my ($m, $s) = stats($x);
  push @out, sprintf "%.3g", $m;
  push @out, sprintf "%.3g", $s;
  my $P = NBBootstrapTest($stat, $x, $ncrit, $maxn);
  #print F "$gene\t$m\t$s\t$P\n";
  printf F "%s\t%.4g\t%.4g\t%.4g\n", $gene, $m, $s, $P;
  Percent(($i-$n1)/($n2-$n1));
}
print "\n";

close F



=head1 OPTIONS

=over 5

=item B<-batch>=I<number>

Batch number.

=item B<-logdir>=I<path>

Path to the directory where all intermediate files are stored.

=item B<-outfile>=I<file>

Output file containing DE results.

=item B<-batchsize>=I<number>

Number of genes in one batch. Default is 30.

=item B<-ncrit>=I<number>

Critical number of positive results to stop the simulation. See function C<NBBootstrapTest> in module C<Distribution.pm>. Default is 30.

=item B<-maxn>=I<number>

Maximum number of simulations. Default is 1e7.

=item B<-stat>=I<m|a>

Statistic to use with Meintanis test. 'm' is for Meintanis, 'a' for Anderson-Darling. The default value is 'm'.


=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut