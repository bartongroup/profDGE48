#!/sw/bin/perl
require Schedule::DRMAAc;
require Cwd;
require PDL;
require PDL::NiceSlice;
require GRNASeq;
require GetOpt;
require Tools;

=head1 NAME

B<grid_launcher_powertest.pl>

=head1 DESCRIPTION

Grid launcher for power tests. This script launches differential expression power tests on the cluster. The DE test of choice (see C<simple_de_test.pl> for details) is run C<nsim> times for increasing number of replicates, starting from 2 and ending at C<maxn>. P-values and fold changes are stored in sqlite databases, to be later processed by C<combine_bs_pvalues.pl> and <make_powerstats_db.pl>.

Note: this script is for simple tests (t-test, log t-test, MW). It does the same job as C<generic_wrapper.py> for tools like edgeR or DEseq.

Requires script C<one_powertest_simple.pl>.

=cut

# Nick: This doesn't actually check if the cluster jobs exit with an exit flag!!
# if they do, you just get an apparent sucess!

use strict;
use warnings;

use Schedule::DRMAAc qw/:all/;
use Cwd;
use PDL;
use PDL::NiceSlice;
use GRNASeq;
use GetOpt;
use Tools;

$| = 1;
Stamp();

my $cwd = getcwd; 

my $test = 't';
my $maxn = 42;
my $nsim = 100;
my $script = "$cwd/scripts/one_powertest_simple.pl";
my %options;
my $logdir = "$cwd/de_tests_logs";
my $dbdir = "$cwd/de_tests_db";
my $queue = '64bit-pri.q,64bit.q,c6145.q';
my $conditions;
my $randcond;
my ($spike, $normspike);
my $defsfile;
my $simpledefile= "$cwd/scripts/simple_de_test.pl";

GetMyOptions(
  'test=s' => \$test,
  'maxn=i' => \$maxn,
  'nsim=i' => \$nsim,
  'options=s' => \%options,
  'conditions=s' => \$conditions,
  'logdir=s' => \$logdir,
  'dbdir=s' => \$dbdir,
  'script=s' => \$script,
  'queue=s' => \$queue,
  'randcond=s' => \$randcond,
  spike => \$spike,
  #normspike => \$normspike,
  #'defsfile=s' => \$defsfile,
  'simpleDEscript=s' => \$simpledefile
);

die "Need valid -script\n" unless defined $script and -e $script;

###############################################

my $Q = join(" ", map {"-q $_"} split /,/, $queue);

my @opt = map {"-$_=$options{$_}"} keys %options;
push @opt, '-spike' if $spike;
#push @opt, '-normspike' if $normspike;
push @opt, "-norm=$norm" if $norm;
@conds = split /,/, $conditions if defined $conditions;

mkdir $logdir, 0777 unless -d $logdir;
mkdir $dbdir, 0777 unless -d $dbdir;

# submit jobs to cluster

my @jobs = ();
my @res = ();
print "Submitting jobs to cluster...       ";

my ($err, $diag) = drmaa_init(undef);
die drmaa_strerror($err) . "\n" . $diag if $err;
$maxn = 20 if defined $randcond;

for my $nrep (2 .. $maxn)
{
  my $job = SubmitBatch($nrep);
  push @jobs, $job;
}
print "\b\b\b\b done    \n";
print "Submitted ", $maxn-1, " jobs to the cluster\n\n";

# monitor cluster activity until all the jobs are finished

my $njobs = @jobs;
my $alldone = 0;
while(!$alldone)
{
  my $queued = 0;
  my $running = 0;
  my $finished = 0;
  my $failed = 0;
  for my $jobid (@jobs)
  {
    my ($err, $ps) = drmaa_job_ps($jobid);
    if($ps == $DRMAA_PS_QUEUED_ACTIVE) {$queued++}
    elsif($ps == $DRMAA_PS_RUNNING) {$running++}
    elsif($ps == $DRMAA_PS_DONE) {$finished++}
    elsif($ps == $DRMAA_PS_FAILED) {$failed++}
  }
  $alldone = ($finished + $failed >= $njobs);
  my $s = sprintf "Queued %3d  Running %3d  Finished %3d  Failed %3d", $queued, $running, $finished, $failed;
  my $b = "\b" x length($s);
  print "$b$s";
  sleep 5 unless $alldone;
}
print "\n";

###############################################################

sub SubmitBatch
{
  my ($nrep) = @_;
  
  my $num = sprintf "%05d", $nrep;

  my $dbfile = (defined $randcond) ? "de_${test}_same${randcond}_rep${nrep}.db" : "de_${test}_rep${nrep}.db";
  my @args = ("-script=$simpledefile", "-nrep=$nrep", "-test=$test", "-nsim=$nsim", "-dbfile=$dbdir/$dbfile", "-defsfile=$defsfile", @opt);
  push @args, "-randcond=$randcond" if defined $randcond;
 
  my ($err, $jt, $diag) = drmaa_allocate_job_template();
  die drmaa_strerror($err) . "\n" . $diag if $err;

  drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $script); 
  drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, "-clear $Q" );  
  drmaa_set_vector_attribute($jt, $DRMAA_V_ARGV, \@args);
  drmaa_set_vector_attribute($jt, $DRMAA_V_ENV, ["PERL5LIB=$ENV{PERL5LIB}", "PATH=$ENV{PATH}"]);
  drmaa_set_attribute($jt, $DRMAA_JOB_NAME, "PT$num");
  drmaa_set_attribute($jt, $DRMAA_OUTPUT_PATH, ':' . $logdir);
  drmaa_set_attribute($jt, $DRMAA_ERROR_PATH, ':' . $logdir);
  drmaa_set_attribute($jt, $DRMAA_WD, $cwd);
  
  ($err, my $jobid, $diag) = drmaa_run_job($jt);
  die drmaa_strerror($err) . "\n" . $diag if $err;
   
  return $jobid
}



=head1 SYNOPSIS

  grid_launcher_powertest.pl -test=st
  grid_launcher_powertest.pl -test=lt -nsim=1000 -logdir=./lt_test -options norm=deseq
  

=head1 OPTIONS

=over 5

=item B<-script>=I<path/file>

File name (with absolute path) containing the script C<one_powertest_simple.pl>. This is a wrapper around the actual DE test script C<simple_de_test.pl>. If not specified, the dafult location, C<./scripts/one_powertest_simple.pl>, will be used.

=item B<-simpleDEscript>=I<path/file>

File name (with absolute path) containing the script C<simple_de_test.pl>. C<one_powertest_simple.pl> is a wrapper around this, the actual DE test script C<simple_de_test.pl>. If not specified, the dafult location, C<./simple_de_test.pl>, will be used.

=item B<-logdir>=I<path>

Path to the directory where all STDOUT and STDERR files from individual cluster jobs are stored. The default value is C<de_tests_logs>.

=item B<-dbdir>=I<path>

Path to the directory where all results are saved. These are Sqlite databases, one file per number of replicates. The default value is C<de_tests_db>.

=item B<-test>=I<[t|lt|st|mw|ks]>

Type of differential expression test.

  t - standard t-test for two samples with identical variances
  lt - logarithmic t-test using approximately normal statistic Z = (log m2 - log m1) / sqrt(((s1/m1)**2)/n1 + ((s2/m2)**2)/n2), where m1, m2 are sample means and s1, s2 are sample standard deviations
  st - shrinkage variance test of Opgen-Rhein and Strimmer (2007)
  mw - Mann-Whitney test
  ks - Kolmogorov-Smirnov test; B<warning> - KS test does not work properly for discrete data!
  
=item B<-maxn>=I<number>

Maximum number of replicates to consider. The scrip will iterate the DE test from 2 replicates up to this number. It should be smaller than the total number of replicates available!

=item B<-nsim>=I<number>

Number of bootstrap runs for the given number of replicates. Median p-values will be calculated over these runs. The default value is 30.

=item B<-spike>

If specified, spike-ins will added to existing genes and spike-in-clean replicates will be selected.

=item B<-norm>=I<string>

Normalization to be used (deseq, totcountm, spikein). If not defined, deseq will be used.

=item B<-defsfile>=I<file>

The full path to the definitions file to be used. Defaults to 'defs.dat'.

=item B<-options>=I<string>

Any additional options passed to the executed script (see Synopsis for usage).

=item B<-queue>=I<string>

Cluster queue name(s), comma-delimited.

=item B<-randcond>=I<string>

If specified, test will be performed on randomly selected replicates from the same condition, defined here.


=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut
