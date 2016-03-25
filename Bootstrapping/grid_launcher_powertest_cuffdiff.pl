#!/sw/bin/perl

=head1 NAME

B<grid_launcher_powertest_cuffdiff.pl>

=head1 DESCRIPTION

This script launches cuffdiff power tests on the cluster. Cuffdiff is run C<nsim> times for increasing number of replicates, starting from 2 and ending at C<maxn>. Fold changes and p-value is calculated for each gene, for the given number of replicates across all simulations. These p-values can be later used to compute power of the test as a function of number of replicates.

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use Schedule::DRMAAc qw/:all/;
use Cwd;
use PDL;
use PDL::NiceSlice;
use GRNASeq;
use GetOpt;
use CompBio::Tools;

$| = 1;
Stamp();

my $created = "Created with $0 " . join(" ", @ARGV) . "\n";

my $cwd = getcwd; 

my $maxn = 42;
my $nsim = 30;
my $totrep = 48;
my $cuffdiff_script = "$cwd/GRNA_Mapping_Tools/run_cuffdiff.pl"; 
my %options;
my $logdir = "$cwd/de_tests_logs";
my $dbdir = "$cwd/de_tests_db";
my ($genlist, $datapath, $genome, $annotation);
my $conditions;
my $replist;
my $queue = 'c6145-long.q';
my $cuffdir = "$cwd/NOBACK/cuffdiff_powertests";
my $spike;
my $threads = 1;
my $email;
my $randcond;
my ($man, $help);

GetOptions(
  'maxn=i' => \$maxn,
  'nsim=i' => \$nsim,
  #'totrep=i' => \$totrep,
  #'options=s' => \%options,
  #'conditions=s' => \$conditions,
  #'dbdir=s' => \$dbdir,
  'datapath=s' => \$datapath,
  'genome=s' => \$genome,
  'annotation=s' => \$annotation,
  'genlist=s' => \$genlist,
  'cuffdiff_script=s' => \$cuffdiff_script,
  #'replist=s' => \$replist,
  'queue=s' => \$queue,
  'cuffdir=s' => \$cuffdir,
  #'threads=i' => \$threads,
  #'email=s' => \$email,
  'randcond=s' => \$randcond,
  #spike => \$spike,
  man => \$man,
  help => \$help
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;

$maxn = 20 if $randcond && $maxn == 42;

die "Need valid -cuffdiff_script\n" unless defined $cuffdiff_script && -e $cuffdiff_script;
die "Need valid -genlist\n" unless defined $genlist && -e $genlist;
die "Need -datapath\n" unless defined $datapath;
die "Need valid -genome\n" unless defined $genome && -e $genome;
die "Need valid -annotation\n" unless defined $annotation && -e $annotation;

###############################################

ReadDefs();

mkdir $cuffdir, 0777 unless -d $cuffdir;
mkdir $dbdir, 0777 unless -d $dbdir;

open R, ">$cuffdir/README" or die "Cannot open log file\n";
print R "$created\n";
close R;

my $m = ($email) ? "-M $email -m a" : '';
my @opt = map {"-$_=$options{$_}"} keys %options;
@conds = split /,/, $conditions if defined $conditions;


my @replist = (defined $replist) ? split(/,/, $replist) : (2 .. $maxn);

# submit jobs to cluster

my @jobs = ();
my @tfiles = ();
my @cmd = ();
print "Submitting jobs to cluster...         ";

my ($err, $diag) = drmaa_init(undef);
die drmaa_strerror($err) . "\n" . $diag if $err;

for my $nrep (@replist)
{
  for my $i (1 .. $nsim)
  {
    my $num = sprintf "r%02d_i%03d", $nrep, $i;
    print "\b\b\b\b\b\b\b\b$num";

    # empirical memory requirement for cuffdiff
    # estimated for yeast, can be much bigger for larger genomes!
    my $mem = ceil(0.32 * $nrep + 1.4);

    my $tfile = "$cuffdir/cuffdiff_$num.tsv";
    my $outdir = "$cuffdir/cuffdiff_out_$num";
    mkdir $outdir, 0777 unless -d $outdir;
    push @tfiles, $tfile;
    my $inrep = (defined $randcond) ? randrepsame($nrep) : randrep($nrep);
    my $job = SubmitBatch($num, $tfile, $outdir, $inrep, $mem);
    push @jobs, $job;
  
    push @cmd, "cuffdiff powertest batch $i, nrep=$nrep, tfile=$tfile, outdir=$outdir, inrep=$inrep";
  }
}
print "\b\b\b\b\b\b\b\b done       \n";
print "Submitted ", scalar @jobs, " jobs to the cluster\n\n";

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

##############################################

sub SubmitBatch
{
  my ($num, $outfile, $outdir, $inrep, $mem) = @_;
  
  my $logdir = $outdir;
  
  my @args = (
    "-datapath=$datapath",
    "-genome=$genome",
    "-annotation=$annotation",
    "-outfile=$outfile",
    "-outpath=$outdir",
    "-inrep=$inrep",
    "-threads=$threads"
  );
  
  push @args, "-conditions=${randcond}A,${randcond}B" if defined $randcond;
  
  my ($err, $jt, $diag) = drmaa_allocate_job_template();
  die drmaa_strerror($err) . "\n" . $diag if $err;

  my $pe = ($threads > 1) ? "-pe smp $threads -R y" : '';

  drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $cuffdiff_script);
  drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, "$m -clear -q $queue -l ram=${mem}G $pe" );  
  #  -pe smp 4 -R y removed 
  drmaa_set_vector_attribute($jt, $DRMAA_V_ARGV, \@args);
  drmaa_set_vector_attribute($jt, $DRMAA_V_ENV, ["PERL5LIB=$ENV{PERL5LIB}", "PATH=$ENV{PATH}"]);
  drmaa_set_attribute($jt, $DRMAA_JOB_NAME, "CD$num");
  drmaa_set_attribute($jt, $DRMAA_OUTPUT_PATH, ':' . $logdir);
  drmaa_set_attribute($jt, $DRMAA_ERROR_PATH, ':' . $logdir);
  drmaa_set_attribute($jt, $DRMAA_WD, $cwd);
  
  ($err, my $jobid, $diag) = drmaa_run_job($jt);
  die drmaa_strerror($err) . "\n" . $diag if $err;
   
  return $jobid
}




sub randrep
{
  my $nrep = shift;
 
  my $ex = ($spike) ? $SpikeinCleanExclude : $CleanExclude; 
  my %e = HashList($ex, 'exclude');
  
  my @sel;
  for my $cond (@conds)
  {
    my $ex = pdl($e{$cond}{exclude}) - 1;
    my $all = sequence($totrep);
    my $sel = setops($all, 'XOR', $ex);
    #print "ex=$ex\nall=$all\nxsel=$sel\n";
    
    my $r = qsorti random $sel;
    $r = $r(0:$nrep - 1);  # random selection
    $sel = qsort $sel($r);
    my @rsel = map {$cond . $_} list($sel + 1);
    my $rsel = join ",", @rsel;
    push @sel, $rsel; 
  }
  my $selected = join ",", @sel;
}


sub randrepsame
{
  my $nrep = shift;
  
  my $exclude = ($spike) ? $SpikeinCleanExclude : $CleanExclude; 
  my %e = HashList($exclude, 'exclude');

  my $ex = pdl($e{$randcond}{exclude}) - 1;
  my $all = sequence($totrep);
  my $xsel = setops($all, 'XOR', $ex);

  my $n = nelem($xsel);
  my $n2 = floor($n / 2);
  die "Number of replicates, $nrep, too large\n" if $nrep > $n2;

  my $r = qsorti random $n;
  my $r1 = $r(0:$nrep-1);  # random selection
  my $r2 = $r($nrep:2*$nrep-1);
  my $xsel1 = $xsel($r1);
  my $xsel2 = $xsel($r2);
  
  my @rsel1 = map {$randcond . 'A' . $_} list($xsel1 + 1);
  my @rsel2 = map {$randcond . 'B' . $_} list($xsel2 + 1);
  
  my $selected = join ",", (@rsel1, @rsel2);
}


=head1 SYNOPSIS

  grid_launcher_powertest_cuffdiff.pl -datapath $g/mapping/genome_biol_reps -annotation Saccharomyces_cerevisiae.EF4.68.gtf -genome Scerevisiae68_ERCC92.fasta -genlist genlist.tsv -nsim=30 -maxn=40 -cuffdiff_script=./run_cuffdiff.pl
    
=head1 OPTIONS

=over 5

=item B<-cuffdiff_script>=I<path/file>

File name (with absolute path) containing the script C<run_cuffdiff.pl>.

=item B<-datapath>=C<path>

Path to the data directory of the experiment. Each sub-directory in this will be treated as a separate condition in the experiment (with the name of the condition matching the name of the directory), and each .bam file in each directory is a replicate for that condition.

=item B<-annotation>=C<filename>

Path to the .gff feature annotation file for the data. This file should match the feature you are counting the RNA-Seq expression for.

=item B<-genome>=C<filename>

Path to the FASTA file with reference genome. This is used by cuffdiff to run its bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates. You should have write access to the directory with this file as cufflinks will attempt to write an index file there. It might be a good idea to copy genome file (or a link) into your local directory.

=item B<-genlist>=C<filename>

Path to the file containing a list of all genes of interest. It can contain multiple columns, the first column should contain (case-insensitive) gene names.

=item B<-maxn>=I<number>

Maximum number of replicates to consider. The scrip will iterate the DE test from 2 replicates up to this number. It should be smaller than the total number of replicates available!

=item B<-nsim>=I<number>

Number of bootstrap runs for the given number of replicates. Median p-values will be calculated over these runs. The default value is 30.

=item B<-cuffdir>=I<path>

Path to a directory where all cuffdiff output files will be stored. It will be divided into subdirectories for each number of replicates and each iteration. These are all intermediate files, so they can be deleted when you are satisfied with the result.

=item B<-randcond>=I<string>

If specified, cuffdiff will be performed on randomly selected replicates from the same condition, defined here.

=item B<-queue>=I<string>

Cluster queue name or names (comma-delimited).

=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut