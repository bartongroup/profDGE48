#!/sw/bin/perl

=head1 NAME

B<grid_launcher_nbtest.pl>

=head1 DESCRIPTION

Runs negative binomial test (script C<one_nb_test.pl>) on the cluster.

=cut

use strict;
use warnings;

use Schedule::DRMAAc qw/:all/;
use Cwd;
use PDL;
use PDL::NiceSlice;
use GRNASeq;
use GetOpt;
use Tools;

ReadDefs();
Stamp();

$| = 1;

my $cwd = getcwd; 

my $script = "$cwd/scripts/one_nb_test.pl";
my $test = 'meintanis';   # mi or meintanis
my $stat = 'm';
my %options;
my $outfile;
my $logdir = "$cwd/NOBACK/nb_test";
my $batchsize = 30;
my $ncrit = 30;
my $maxn = 1e7;
my $queue = 'c6100.q,64bit-pri.q';
my $mem = '4000M';
my $maxmem = '32000M';
GetMyOptions(
  'options=s' => \%options,
  'outfile=s' => \$outfile,
  #'test=s' => \$test,
  'stat=s' => \$stat, 
  'logdir=s' => \$logdir,
  'script=s' => \$script,
  'batchsize=s' => \$batchsize,
  'ncrit=i' => \$ncrit,
  'maxn=f' => \$maxn,
  'queue=s' => \$queue,
  'mem=s' => \$mem
);



die "Need -script\n" unless defined $script and -e $script;
die "Need -outfile\n" unless defined $outfile;

###############################################

my $cl = '';
$cl = '_clean' if $clean;
$cl = '_spclean' if $spclean;
my $root = "${cond}_${type}_${norm}${cl}_${stat}";

my @opt = map {"-$_=$options{$_}"} keys %options;
mkdir $logdir, 0777 unless -d $logdir;

push @opt, ("-cond=$cond", "-norm=$norm", "-nonzero=$nonzero", "-type=$type", "-test=$test");
push @opt, "-clean" if $clean;
push @opt, "-ncrit=$ncrit" if $ncrit;

push @opt, "-maxn=$maxn" if $maxn;


my $ngen = ngenes();
my $nbatch = ceil($ngen / $batchsize);

my @jobs = ();
my @res = ();
print "Submitting jobs to cluster...      ";

my ($err, $diag) = drmaa_init(undef);
die drmaa_strerror($err) . "\n" . $diag if $err;

for my $batch (1 .. $nbatch)
{
  my $num = sprintf "%05d", $batch;
  my $resfile = "$logdir/${root}_$num.txt";
  printf "\b\b\b\b\b$num";

  my @args = ["-batch=$batch", "-batchsize=$batchsize", "-outfile=$resfile", "-stat=$stat", @opt];
 
  my ($err, $jt, $diag) = drmaa_allocate_job_template();
  die drmaa_strerror($err) . "\n" . $diag if $err;

  drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, "-clear -q $queue -l ram=$mem -l h_vmem=$maxmem" );
  drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $script); 
  drmaa_set_vector_attribute($jt, $DRMAA_V_ARGV, @args);
  drmaa_set_attribute($jt, $DRMAA_JOB_NAME, "NB$num");
  drmaa_set_attribute($jt, $DRMAA_OUTPUT_PATH, ':' . $logdir);
  drmaa_set_attribute($jt, $DRMAA_ERROR_PATH, ':' . $logdir);
  drmaa_set_vector_attribute($jt, $DRMAA_V_ENV, ["PERL5LIB=$ENV{PERL5LIB}", "PATH=$ENV{PATH}"]);
  drmaa_set_attribute($jt, $DRMAA_WD, $cwd);
  
  ($err, my $jobid, $diag) = drmaa_run_job($jt);
  die drmaa_strerror($err) . "\n" . $diag if $err;
  
  push @jobs, $jobid;
  push @res, $resfile;
}

print "\b\b\b\b\b done    \n";
print "Submitted $nbatch jobs with $ngen genes to the cluster\n\n";

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
  sleep 10 unless $alldone;
}
print "\n";

# Aggregate results

print "Aggregating results...";
open OUT, ">$outfile" or die "Cannot open output file $outfile\n";
for my $r (@res)
{
  open F, $r or die "Cannot open $r\n";
  while(<F>)
    {print OUT}
}
print " done\n";

# Cleaning up

for my $file (@res)
  {unlink $file}

print "Results in $outfile\n";








sub ngenes
{
  my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel) = GetNormalizedData(
    condition=>$cond,
    replicate=>$rep,
    exclude=>$exclude,
    nonzero=>1,
    norm=>$norm,
    type=>$type,
    clean=>$clean,
    quiet=>1
  );
  return $ngen;
}


###############################################################




=head1 SYNOPSIS

  grid_launcher_nbtest.pl -script=./one_nb_test.pl -logdir=./nb_test -outfile=nbtest.txt

=head1 OPTIONS

=over 5

=item B<-script>=I<path/file>

Location (full path and name) of the script C<one_nb_test.pl>.

=item B<-logdir>=I<path>

Path to the directory where all intermediate files are stored.

=item B<-outfile>=I<file>

Aggregated output file containing DE results.

=item B<-batchsize>=I<number>

Number of genes in one batch. Default is 30.

=item B<-ncrit>=I<number>

Critical number of positive results to stop the simulation. See function C<NBBootstrapTest> in module C<Distribution.pm>. Default is 30.

=item B<-maxn>=I<number>

Maximum number of simulations. Default is 1e7.


=item B<-stat>=I<m|a>

Statistic to use with Meintanis test. 'm' is for Meintanis, 'a' for Anderson-Darling. The default value is 'm'.

=item B<-options>=I<string>

Any other options passed to the script.

=item B<-queue>=I<string>

Cluster queue name(s), comma delimited.

=item B<-mem>=I<value>

Memory requirement for each job. The default value is '4000M'.

=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut