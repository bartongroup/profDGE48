#!/sw/bin/perl

=head1 NAME

B<grid_launcher_DE.pl>

=head1 DESCRIPTION

This is a wrapper script that divides gene counts into small chunks (default 30 genes each) and runs DE scripts on the cluster, using DRMAA. The results are then aggregated into one file. This script is required to calculate permutation and bootstrap tests.

=cut

use strict;
use warnings;

use Schedule::DRMAAc qw/:all/;
use Cwd;
use PDL;
use PDL::NiceSlice;
use GRNASeq;
use GetOpt;
use CompBio::Tools;

$| = 1;
Stamp();

my @conds = qw(WT Snf2);
my $cwd = getcwd; 

my $script;
my %options;
my $outfile;
my $logdir = "$cwd/grid_output";
my $conditions;
my $batchsize = 30;
GetMyOptions(
  'options=s' => \%options,
  'outfile=s' => \$outfile,
  #'conditions=s' => \$conditions,
  'logdir=s' => \$logdir,
  'script=s' => \$script,
  'batchsize=s' => \$batchsize
);

die "Need -script\n" unless defined $script and -e $script;
die "Need -outfile\n" unless defined $outfile;

###############################################

my @opt = map {"-$_=$options{$_}"} keys %options;
mkdir $logdir, 0777 unless -d $logdir;
@conds = split /,/, $conditions if defined $conditions;

# Read gene list from the first condition file

my $genlist = CountFile('yeast', $conds[0], undef, 'raw');
my @genes = ReadGeneList($genlist);

# Read count data

my %data = ();
my %g2i = ();
for my $condition (@conds)
{
  my @dat = GetNormalizedData(
    condition=>$condition,
    exclude=>$exclude,
    include=>$include,
    nonzero=>$nonzero,
    norm=>$norm,
    type=>$type,
    clean=>$clean,
    spclean=>$spclean
  );
  $data{$condition} = \@dat;
  my %g = GenHash($dat[5], $dat[7]);
  $g2i{$condition} = \%g;
  # @dat is  ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel) 
}

# submit jobs to cluster

my $cond1 = $conds[0];
my $cond2 = $conds[1];
my $d1 = $data{$cond1}[0];
my $d2 = $data{$cond2}[0];

my $buffer = '';
my @jobs = ();
my @res = ();
my @bat = ();
my $ngen = 0;
print "Submitting jobs to cluster...       ";

my ($err, $diag) = drmaa_init(undef);
die drmaa_strerror($err) . "\n" . $diag if $err;

my $batch = 1;
my $n = 0;
for my $gene (@genes)
{
  $gene =~ tr/A-Z/a-z/;
  my $i1 = $g2i{$cond1}{$gene};
  my $i2 = $g2i{$cond2}{$gene};
  if(defined $i1 && defined $i2)
  {
    $ngen++;
    
    my $x1 = $d1(,$i1;-);
    my $x2 = $d2(,$i2;-);
    
    my ($m1) = stats($x1);
    my ($m2) = stats($x2);
    my $fc = $m2 / $m1;
    
    my $s1 = join "\t", list($x1);
    my $s2 = join "\t", list($x2);
    
    $buffer .= "$gene\n$s1\n$s2\n";
    $n++;
    
    if($n >= $batchsize)
    {
      printf "\b\b\b\b%4d", $batch;
      my ($job, $batchfile, $resfile) = SubmitBatch($batch, $buffer);
      push @jobs, $job;
      push @res, $resfile;
      push @bat, $batchfile;
      $buffer = '';
      $n = 0;
      $batch++;
    }
  }
  #last if $batch > 10;
}
unless($buffer eq '')  # last one
{
  my ($job, $batchfile, $resfile) = SubmitBatch($batch, $buffer);
  push @jobs, $job;
  push @res, $resfile;
  push @bat, $batchfile;
  $batch++
}

$batch--;
print "\b\b\b\b done    \n";
print "Submitted $batch jobs with $ngen genes to the cluster\n\n";

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

for my $file (@res, @bat)
  {unlink $file}

print "Results in $outfile\n";

###############################################################

sub SubmitBatch
{
  my ($batchnum, $str) = @_;
  
  my $num = sprintf "%05d", $batchnum;
  my $batchfile = "$logdir/batch_${$}_$num.txt";
  my $resfile = "$logdir/results_${$}_$num.txt";
  local *F;
  open F, ">$batchfile" or die "$!\n";
  print F $str;
  close F;

  my @args = ["-infile=$batchfile", "-outfile=$resfile", @opt];
 
  my ($err, $jt, $diag) = drmaa_allocate_job_template();
  die drmaa_strerror($err) . "\n" . $diag if $err;

  drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $script); 
  drmaa_set_vector_attribute($jt, $DRMAA_V_ARGV, @args);
  drmaa_set_vector_attribute($jt, $DRMAA_V_ENV, ["PERL5LIB=$ENV{PERL5LIB}", "PATH=$ENV{PATH}"]);
  drmaa_set_attribute($jt, $DRMAA_JOB_NAME, "DE$num");
  drmaa_set_attribute($jt, $DRMAA_OUTPUT_PATH, ':' . $logdir);
  drmaa_set_attribute($jt, $DRMAA_ERROR_PATH, ':' . $logdir);
  drmaa_set_attribute($jt, $DRMAA_WD, $cwd);
  
  ($err, my $jobid, $diag) = drmaa_run_job($jt);
  die drmaa_strerror($err) . "\n" . $diag if $err;
   
  return $jobid, $batchfile, $resfile
}



=head1 SYNOPSIS

  grid_launcher_DE.pl -script=./one_bootstrap_ET93.pl -logdir=./boot_output -outfile=boot_result.tsv -options nboot=100000
  grid_launcher_DE.pl -script=./one_permutest.pl -logdir=./perm_output -outfile=perm_result.tsv
  
=head1 OPTIONS

=over 5

=item B<-script>=I<path/file>

File name (with absolute path) containing the script to run DE. The script should have two mandatory options, C<-infile> and C<-outfile>, and any other additional options.

Input file is a chunk of gene counts in individual replicates and both conditions. There are three lines per each gene: gene name, counts for condition one (tab-delimited) and counts for conditions two. Such file is created by this script and passed on to the DE script.

Output file countains resulting DE calles. There should be at least five tab-delimited columns in it: gene name, log2 fold change, p-value, mean of condition 1, mean of condition 2. Individual output files are aggregated into the final output file.

Any additional parameters can be included. These are passed on by C<-options> option.

DE scripts available at this moment are C<one_bootstrap_ET93.pl> and C<one_permutest.pl>.

=item B<-logdir>=I<path>

Path to the directory where all intermediate files are stored.

=item B<-outfile>=I<file>

Aggregated output file containing DE results.

=item B<-batchsize>=I<number>

Number of genes in one batch. Default is 30.

=item B<-options>=I<options>

All other options required by the DE script, for example, number if iterations. See synopsis for an example. 

=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut