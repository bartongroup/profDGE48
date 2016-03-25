#!/sw/bin/perl

=head1 NAME

B<make_powerstats_db.pl>

=head1 DESCRIPTION

Calculates power statistics (true positives, true negatives, false positives, false negatives) for each test/tool against a chosen "true standard". The power stats are calculated for a range of replicates, while the true standard is a test result for the full clean data set. You need to run all power bootstraps first. The results are stored in a directory to be used later to create comparative power plots.

This script uses data directly from db files, and counts true/false positives/negatives for each bootstrap separately, then combining these results and reporting mean and standard deviation. An alternative approach is used in C<make_powerstats.pl>, where median p-values across bootstraps are used.

This script is a wrapper around C<one_bs_powerstats.pl>, and it runs it on the cluster.

=cut


use strict;
use warnings;

use DBI;
use Schedule::DRMAAc qw/:all/;
use File::Temp qw(tempfile mktemp);
use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;

use GRNASeq;
use Cwd;

use Tools;
use Stats;

Stamp();
my $cwd = getcwd; 
$| = 1;

my $script = "DE_tool_comparison/one_bs_powerstats.pl";

my $test = 't';
my $truetests = 'mw';
#my $truetests = '';
#my $dbdir = "$cwd/de_tests_db";
my $logdir = "$cwd/powerstats_logs";
#my $dbform = 'de_t_rep%d.db';
my $multicor = 'bh';
my $testinfofile = './de_tests.txt';
my $powerdir = './powerstats_db_ref';
my $maxn = 40;
my $alpha = 0.05;
my $queue = 'c6100.q,64bit-pri.q,64bit.q,c6145.q';
my $genlist = 'genlist.tsv';
my $reffcfile = 'WT_Snf2_deseq_clean_mw_test.tsv';
my $colp;
my ($help, $man);

GetOptions(
  'test=s' => \$test,
#  'dbdir=s' => \$dbdir,
#  'dbform=s' => \$dbform,
  'testinfofile=s' => \$testinfofile,
  'logdir=s' => \$logdir,
  'script=s' => \$script,
  #'alpha=f' => \$alpha,
  'maxn=i' => \$maxn,
  'queue=s' => \$queue,
  'powerdir=s' => \$powerdir,
  'genlist=s' => \$genlist,
  'reffcfile=s' => \$reffcfile,
  #'colp=i' => \$colp,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;


mkdir $logdir, 0777 unless -d $logdir;
mkdir $powerdir, 0777 unless -d $powerdir;


my $Q = join(" ", map {"-q $_"} split /,/, $queue);

###############################################

print "Submitting jobs to cluster...      ";
my @jobs = ();
my ($err, $diag) = drmaa_init(undef);
die drmaa_strerror($err) . "\n" . $diag if $err;

my $nbatch = 0;

my @truetests = split(/\,/, $truetests);
push @truetests, 'self';

my %tests = ReadTestfileInfo($testinfofile);
die "Unrecognized test $test\n" unless defined $tests{$test};
die "No repfiles for test $test\n" if $tests{$test}{repfile} eq '-';

my %jobdat = ();
for my $truetest (@truetests)
{ 
  $truetest = $test if $truetest eq 'self';
  die "Unrecognized test $truetest\n" unless defined $tests{$truetest};

  my %t = %{$tests{$test}};
  my %tt = %{$tests{$truetest}};

  my ($tnsig, $tnnonsig) = FullNSig(%t);

  my @reps = ();
  for my $nrep (2 .. $maxn)
  {
    #my $dbfile = sprintf "$dbdir/$dbform", $nrep;
    my $dbfile = sprintf $t{repfile}, $nrep;
    $dbfile =~ s/tsv$/db/;     # hm...
    next unless -e $dbfile;
  
    $nbatch++;
    printf "\b\b\b%3d", $nbatch;
    my $totfile = "$logdir/" . mktemp("pstat_tot_${test}_${truetest}_n${nrep}_XXXXXX");
    my $fcfile1 = "$logdir/" . mktemp("pstat_fc1_${test}_${truetest}_n${nrep}_XXXXXX");
    my $fcfile2 = "$logdir/" . mktemp("pstat_fc2_${test}_${truetest}_n${nrep}_XXXXXX");
    my $rocfile = "$logdir/" . mktemp("pstat_roc_${test}_${truetest}_n${nrep}_XXXXXX");
    my $nsigfile = "$logdir/" . mktemp("pstat_nsig_${test}_${truetest}_n${nrep}_XXXXXX");
    my $propfile = "$logdir/" . mktemp("pstat_prop_${test}_${truetest}_n${nrep}_XXXXXX");
    push @reps, [$nrep, $totfile, $fcfile1, $fcfile2, $rocfile, $nsigfile, $propfile];

    my $jobname = sprintf "PS%03dR%02d", $nbatch, $nrep;
    my @args = (
      "-dbfile=$dbfile",
      "-totfile=$totfile",
      "-fcfile1=$fcfile1",
      "-fcfile2=$fcfile2",
      "-rocfile=$rocfile",
      "-nsigfile=$nsigfile",
      "-sigpropfile=$propfile",
      "-truefile=$tt{fullfile}",
      "-truecor=$tt{adjustment}",
      "-dbcor=$t{adjustment}",
      "-fcsig=$t{FCsig}",
      "-genlist=$genlist"
    );
    push @args, "-reffcfile=$reffcfile" if defined $reffcfile;
    push @args, "-colp=$colp" if defined $colp;
	#print $script," ",join(" ",@args),"\n";
	 
    my ($err, $jt, $diag) = drmaa_allocate_job_template();
    die drmaa_strerror($err) . "\n" . $diag if $err;
  
    drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, "-clear $Q" );
    drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $script); 
    drmaa_set_vector_attribute($jt, $DRMAA_V_ARGV, \@args);
    drmaa_set_attribute($jt, $DRMAA_JOB_NAME, $jobname);
    drmaa_set_attribute($jt, $DRMAA_OUTPUT_PATH, ':' . $logdir);
    drmaa_set_attribute($jt, $DRMAA_ERROR_PATH, ':' . $logdir);
    drmaa_set_vector_attribute($jt, $DRMAA_V_ENV, ["PERL5LIB=$ENV{PERL5LIB}"]);
    drmaa_set_attribute($jt, $DRMAA_WD, $cwd);
    
    ($err, my $jobid, $diag) = drmaa_run_job($jt);
    die drmaa_strerror($err) . "\n" . $diag if $err;
    
    push @jobs, $jobid;
  }
  
  $jobdat{$truetest} = [\@reps, $tnsig, $tnnonsig];
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
  sleep 1 unless $alldone;
}
print "\n";

print "Aggregating results\n";
for my $truetest (@truetests)
{ 
  my $root = "$powerdir/power_${test}_true-${truetest}";
  my $totfile = "$root.stats";  
  my $fcfile1 = "${root}_fc1.stats";
  my $fcfile2 = "${root}_fc2.stats";
  my $rocfile = "${root}_roc.stats";
  my $nfile = "${root}_nsig.stats";
  my $propfile = "${root}_sigprop.stats";
  my ($reps, $nsig, $nnonsig) = @{$jobdat{$truetest}};
  my @reps = @$reps;
  
  open PROP, ">$propfile" or die;
  open NSIG, ">$nfile" or die;
  open TOT, ">$totfile" or die;
  open FC1, ">$fcfile1" or die;
  open FC2, ">$fcfile2" or die;
  open ROC, ">$rocfile" or die;
  print PROP "#nrep\tGene\tFC\tSignificant fraction\n";
  print NSIG "#nrep\tNsig\n";
  print TOT "#nrep\tTP\tTN\tFP\tFN\tsTP\tsTN\tsFP\tsFN\n";
  print FC1 "#nrep\tfclim\tTP\tTN\tFP\tFN\tsTP\tsTN\tsFP\tsFN\n";
  print FC2 "#nrep\tfclim\tTP\tTN\tFP\tFN\tsTP\tsTN\tsFP\tsFN\n";
  print ROC "#nrep\tfclim\tFPR\tTPR\tsFPR\tsTPR\n";
  for my $r (@reps)
  {
    my ($nrep, $tot, $fc1, $fc2, $roc, $ns, $pr) = @$r;
    unless(-e $tot && -e $fc1 && -e $fc2 && -e$roc && -e $ns)
    {
       print "TMP files missing for n=$nrep\n";
       next
    }
    
    open N, $ns or die "Cannot open tmp file $ns\n";
    my $line = <N>;
    print NSIG "$nrep\t$line";
    close N;
    unlink $ns;
    
    open P, $pr or die "Cannot open tmp file $pr\n";
    while(my $line = <P>)
      {print PROP "$nrep\t$line"}
    close P;
    unlink $pr;
    
    
    open T, $tot or die "Cannot open tmp file $tot\n";
    $line = <T>;
    print TOT "$nrep\t$line";
    close T;
    unlink $tot;
    
    open F1, $fc1 or die "Cannot open tmp file $fc1\n";
    while(my $line = <F1>)
      {print FC1 "$nrep\t$line"}
    close F1;
    unlink $fc1;
    
    open F2, $fc2 or die "Cannot open tmp file $fc2\n";
    while(my $line = <F2>)
      {print FC2 "$nrep\t$line"}
    close F2;
    unlink $fc2;

    open R, $roc or die "Cannot open tmp file $roc\n";
    while(my $line = <R>)
      {print ROC "$nrep\t$line"}
    close R;
    unlink $roc;
    
  }
  close ROC;
  close FC1;
  close FC2;
  close TOT;
  close PROP;
  print "Created files wih root $root\n";
}
print " done\n";  
  


sub FullNSig
{
  my (%test) = @_;
  
  my $file = $test{fullfile};
  die "Cannot find file $file\n" unless defined $file && -e $file;
  my $cor = $test{adjustment};
  
  my $p = rcols $file, 2;
  
  # find total number of "real" positives and negatives
  my $mlim = MulticorLimit($p, $cor, $multicor, $alpha);
  my $nsig = nelem($p->where($p <= $mlim));
  my $nnonsig = nelem($p) - $nsig;
  
  #print $p(0:100), "\n";
  #print "$test{name}  $cor  $multicor  $mlim  $nsig\n";die;
  
  return $nsig, $nnonsig
}

=head1 SYNOPSIS

  make_powerstats_db.pl -test=edger

=head1 OPTIONS

=over 5

=item B<-test>=I<string>

Name of the test/tool to plot results for. As defined in the metadata file.

=item B<-testinfofile>=I<pathname>

A file with DE tools metadata. See C<plot_fcp.pl> for details. Default name is C<de_tests.txt>

=item B<-script>=I<path/file>

This should be the full path to the C<one_bs_powerstats_db.pl> script that this wrapper runs.

=item B<-logdir>=I<path>

Path to the directory where all intermediate files are stored. Default is ./powerstats_logs.

=item B<-powerdir>=I<path>

Path to the directory where results are saved. Default value is powerstats_db_ref.

=item B<-maxn>=I<number>

Maximum number of replicates to be used. The default value is 40.

=item B<-queue>=I<string>

Grid engine queue name (or names, comma-separated).

=item B<-genlist>=I<pathname>

A file with the list of genes. Created by C<make_gene_list.pl>. The default value is C<genlist.tsv>.

=item B<-reffcfile>=I<pathname>

A file with reference log2-fold-changes (in the second colum) and gene names (first column). These fold changes will be used for consistency. This file can be created, e.g., by the following command.

   simple_de_test.pl -test=mw -clean -outfile=WT_Snf2_deseq_clean_mw_test.tsv

=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut