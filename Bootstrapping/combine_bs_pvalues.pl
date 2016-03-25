#!/sw/bin/perl

=head1 NAME

B<combine_bs_pvalues.pl>

=head1 DESCRIPTION

Calculates median of p-values across bootstrap runs. The input is all sqlite .db files created by DE bootstrap (for a range of replicate numbers). The script will create a single .tsv file per replicate number, containing median p-values.

This script is a wrapper around C<one_bs_pvalue.pl>, and it runs it on the cluster.

=cut


use strict;
use warnings;

use DBI;
use Schedule::DRMAAc qw/:all/;

use PDL;
use PDL::NiceSlice;

use GRNASeq;
use GetOpt;
use Cwd;

use CompBio::Tools;
use CompBio::Statistics;

my $cwd = getcwd; 
$| = 1;
Stamp();

my $script = "$cwd/scripts/one_bs_pvalue.pl";

my $dbdir = "$cwd/de_tests_db";
my $logdir = "$cwd/de_tests_logs";
my $dbform = 'de_t_rep%d.db';
my $outform;
my $maxn = 48;
my $logit;
my $queue = '64bit-pri.q,64bit.q,c6145.q';
GetMyOptions(
  'dbdir=s' => \$dbdir,
  'dbform=s' => \$dbform,
  'outform=s' => \$outform,
  'logdir=s' => \$logdir,
  'script=s' => \$script,
  #'maxn=i' => \$maxn,
  'queue=s' => \$queue,
# logit => \$logit          # samR reported FC instead of log2(FC): no longer needed
);

mkdir $logdir, 0777 unless -d $logdir;

unless(defined $outform)
{
  $outform = $dbform;
  $outform =~ s/db$/tsv/;
}

my $Q = join(" ", map {"-q $_"} split /,/, $queue);

###############################################

print "Submitting jobs to cluster...      ";
my @jobs = ();
my ($err, $diag) = drmaa_init(undef);
die drmaa_strerror($err) . "\n" . $diag if $err;

my $nbatch = 0;
for my $nrep (2 .. $maxn)
{
  printf "\b\b%2d", $nrep;
  my $dbfile = sprintf "$dbdir/$dbform", $nrep;
  next unless -e $dbfile;
  
  $nbatch++;
  my $outfile = sprintf "$dbdir/$outform", $nrep;

  my $num = sprintf "%02d", $nrep;
  my @args = ("-dbfile=$dbfile", "-outfile=$outfile");
  push @args, "-logit" if $logit;
 
  my ($err, $jt, $diag) = drmaa_allocate_job_template();
  die drmaa_strerror($err) . "\n" . $diag if $err;

  drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, "-clear $Q" );
  drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $script); 
  drmaa_set_vector_attribute($jt, $DRMAA_V_ARGV, \@args);
  drmaa_set_attribute($jt, $DRMAA_JOB_NAME, "MBS$num");
  drmaa_set_attribute($jt, $DRMAA_OUTPUT_PATH, ':' . $logdir);
  drmaa_set_attribute($jt, $DRMAA_ERROR_PATH, ':' . $logdir);
  drmaa_set_vector_attribute($jt, $DRMAA_V_ENV, ["PERL5LIB=$ENV{PERL5LIB}"]);
  drmaa_set_attribute($jt, $DRMAA_WD, $cwd);
  
  ($err, my $jobid, $diag) = drmaa_run_job($jt);
  die drmaa_strerror($err) . "\n" . $diag if $err;
  
  push @jobs, $jobid;
}

print "\b\b\b\b\b done    \n";
print "Submitted $nbatch jobs to the cluster\n\n";

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
  sleep 3 unless $alldone;
}
print "\n";



=head1 SYNOPSIS

  combine_bs_pvalues.pl -dbdir=./bs_results -dbform=de_ttest_rep%d.db

=head1 OPTIONS

=over 5

=item B<-dbdir>=I<path>

Path to the directory containing .db files.

=item B<-dbform>=I<string>

Format of .db file name, e.g. de_test_rep%d.db, for files de_test_rep2.db, de_test_rep3.db, ..., de_test_rep30.db. The replicates will be cycled between 2 and 48, non-existing files ignored.

=item B<-script>=I<path/file>

This should be the full path to the C<one_mean_bs_pvalue.pl> script that this wrapper runs.

=item B<-logdir>=I<path>

Path to the directory where all intermediate files are stored. Default is ./de_tests_logs.

=item B<-outform>=I<string>

Format for output .tsv files.

=item B<-queue>=I<string>

Cluster queue name or names (comma-separated).

=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut