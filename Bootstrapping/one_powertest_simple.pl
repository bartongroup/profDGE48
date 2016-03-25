#!/sw/bin/perl

=head1 NAME

B<one_powertest_simple.pl>

=head1 DESCRIPTION

Script to perform DE power test, used by C<grid_launcher_powertest.pl>. This is essentially a wrapper around C<simple_de_test.pl>. This script is not supposed to be ran on its own.

=cut


use strict;
use warnings;

use DBI;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;
use File::Temp qw(tempfile mktemp);
use Getopt::Long;

use Time::HiRes qw(gettimeofday tv_interval);

use GRNASeq;
use BootstrapDB;
use Tools;
use Pod::Usage;


$| = 1;


my $created = "Created with $0 " . join(" ", @ARGV) . "\n";

my $script = "simple_de_test.pl";
my ($nrep, $test, $nsim, $dbfile, $spike, $normspike, $norm, $randcond);
my $defsfile;
my ($help, $man);
GetOptions(
  'script=s' => \$script,
  'nrep=i' => \$nrep,
  'randcond=s' => \$randcond,
  'test=s' => \$test,
  'nsim=i' => \$nsim,
  'dbfile=s' => \$dbfile,
  #spike => \$spike,
  #normspike => \$normspike,
  'norm=s' => \$norm,
  #'defsfile=s' => \$defsfile,
  help => \$help,
  man => \$man,
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;


unlink $dbfile if -e $dbfile;

ReadDefs($defsfile);

###############################################

my $type = ($spike) ? 'all' : 'raw';

my $genlist = CountFile('yeast', $conds[0], undef, $type);
my @genes = ReadGeneList($genlist);

my $topdir = (defined $ENV{TMPDIR}) ? $ENV{TMPDIR} : '/tmp';
my $dir = "$topdir/powertests";
mkdir $dir, 0777 unless -d $dir; 

print "Running test $test for $nrep replicates\n";
print "Using temp directory $dir\n";
my @tmp = ();
my @cmd = ();
my @time = ();
for my $i (1 .. $nsim)
{
  my $tmpfile = mktemp("$dir/powertestXXXXXX"); 
  push @tmp, $tmpfile;
  
  #my $c = ($spike) ? '-type=all -spclean -norm=spikein' : '-clean';
  my $c = ($spike) ? '-type=all -spclean' : '-clean';
  #my $c = ($spike) ? '-type=all -clean' : '-clean';
  #my $n = ($normspike) ? '-norm=spikein' : '';
  my $n = "-norm=$norm" if $norm;
  my $r = (defined $randcond) ? "-randrep2=$nrep -randcond=$randcond" : "-randrep=$nrep";
  my $cmd = "$script $c $n $r -test=$test -outfile=$tmpfile -defsfile=$defsfile";
  push @cmd, $cmd;
  my $t0 = [gettimeofday];
  `$cmd`;
  my $elapsed = tv_interval ( $t0, [gettimeofday]);
  push @time, $elapsed;
  Percent($i/$nsim);
}
print "\n";


print "Creating database $dbfile...";
my $db = bsConnect($dbfile);
bsCreateTables($db);
bsFillFeatures($db, \@genes);
print " done\n";

print "Aggregating results...";
bsAggregate($db, \@tmp, \@cmd, \@time, $created);
print " done\n";

=head1 OPTIONS

=over 4

=item B<-test>=I<[t|lt|st|mw|ks]>

Type of differential expression test.

  t - standard t-test for two samples with identical variances
  lt - log-ratio t-test using approximately normal statistic Z = (log m2 - log m1) / sqrt(((s1/m1)**2)/n1 + ((s2/m2)**2)/n2), where m1, m2 are sample means and s1, s2 are sample standard deviations
  st - shrinkage variance test of Opgen-Rhein and Strimmer (2007)
  mw - Mann-Whitney test
  ks - Kolmogorov-Smirnov test; B<warning> - KS test does not work properly for discrete data!
  
 =item B<-nrep>=I<number>
 
 Number of replicates to run the test on.
  

=item B<-randcond>=I<string>

If specified, test will be performed on randomly selected replicates from the same condition, defined here.

=item B<-nsim>=I<number>
 
 Number of bootstraps.
 
=item B<-dbfile>=I<pathname>

Name of the sqlite .db file to create. 


=item B<-norm>=I<string>

Expression data normalization. Input expression files contain raw counts. These are normalized on the fly, using one of the following methods:

  none - no normalization
  deseq - DESeq normalization (default)
  tmm - trimmed mean of M values normalization
  fpkm - approximate FPKM normalization (not identical to cuffdiff!)
  totcount - total count
  totcountm - total count with normalization factors scaled to the mean of 1



=back

