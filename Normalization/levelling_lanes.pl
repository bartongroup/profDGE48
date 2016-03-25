#!/sw/bin/perl

=head1 NAME

B<levelling_lanes.pl>  

Equal count (levelling) normalization of fastq files for lanes. Uses the grid engine.

=cut


use strict;
use warnings;

use Schedule::DRMAAc qw/:all/;
use PDL;
use PDL::NiceSlice;
use Cwd;

use CompBio::Tools;
use Getopt::Long;
use Pod::Usage;

use GRNASeq;
our $topdir;


$| = 1;
my $cwd = getcwd;

my $levdir = "$topdir/raw/levelled_lanes";
my @conds = qw(WT Snf2);
my $logdir = "$cwd/NOBACK/levelling";
my $script = "$cwd/local_scripts/one_levelling.pl";
my $queue = '64bit-pri.q,64bit.q,c6145.q';

my ($help, $man);
GetOptions(
  'levdir=s' => \$levdir,
  'logdir=s' => \$logdir,
  'script=s' => \$script,
  'queue=s' => \$queue,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;



my @run = ();
my ($cnt, $mins) = ReadCountFiles();
for my $cond (keys %$cnt)
{
  for my $rep (sort {$a <=> $b} keys %{$cnt->{$cond}})
  {
    my $min = $mins->{$cond}{$rep};
    for my $lane (sort {$a <=> $b} keys %{$cnt->{$cond}{$rep}})
    {
      my ($infile, $count) = @{$cnt->{$cond}{$rep}{$lane}};
      my $outfile;
      if($infile =~ /raw_reads\/(\d+)\/(.*)\.sanfastq\.gz/)
      {
        my $num = $1;
        my $root = $2;
        my $dir = "$levdir/$num";
        mkdir $dir, 0777, unless -d $dir;
        $outfile = "$dir/${root}_levelled.fastq.gz"
      }
      else
        {die "Unrecognized file name $infile\n"}
      #print "$cond $rep $lane $count  ";
      #print "$infile => $outfile\n";
      #Level($infile, $outfile, $count, $min);
      push @run, [$infile, $outfile, $count, $min];
      #print "\n";
    }
  }
}



my ($err, $diag) = drmaa_init(undef);
die drmaa_strerror($err) . "\n" . $diag if $err;

print "Submitting jobs to the cluster...";
my $n = 0;
for my $r (@run)
{
  my ($infile, $outfile, $count, $min) = @$r;
  
  my $num = sprintf "%03d", $n;
  my $resfile = "$logdir/LEV$num.txt";
  #printf "\b\b\b\b\b$num";

  my @args = ["-infile=$infile", "-outfile=$outfile", "-count=$count", "-min=$min"];
 
  my ($err, $jt, $diag) = drmaa_allocate_job_template();
  die drmaa_strerror($err) . "\n" . $diag if $err;

  drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, "-clear -q $queue" );
  drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $script); 
  drmaa_set_vector_attribute($jt, $DRMAA_V_ARGV, @args);
  drmaa_set_attribute($jt, $DRMAA_JOB_NAME, "NB$num");
  drmaa_set_attribute($jt, $DRMAA_OUTPUT_PATH, ':' . $logdir);
  drmaa_set_attribute($jt, $DRMAA_ERROR_PATH, ':' . $logdir);
  drmaa_set_vector_attribute($jt, $DRMAA_V_ENV, ["PERL5LIB=$ENV{PERL5LIB}", "PATH=$ENV{PATH}"]);
  drmaa_set_attribute($jt, $DRMAA_WD, $cwd);
  
  ($err, my $jobid, $diag) = drmaa_run_job($jt);
  die drmaa_strerror($err) . "\n" . $diag if $err;
  
  $n++
}
print " done\nSubmitted $n jobs.\n";



##############################################

sub ReadCountFiles
{
  local *F;
  my $cfile = "all_readcount_lane.dat";
  
  my %cnt = ();
  open F, $cfile or die;
  while(my $line = <F>)
  {
    chomp $line;
    my ($cond, $rep, $lane, $file, $count) = split " ", $line;
    $cnt{$cond}{$rep}{$lane} = [$file, $count];
  }
  
  my %mins = ();
  for my $cond (keys %cnt)
  {
    my %x = %{$cnt{$cond}};
    for my $rep (keys %x)
    {
      my %y = %{$x{$rep}};
      my $min = 1e32;
      for my $lane (keys %y)
        {$min = $y{$lane}[1] if $y{$lane}[1] < $min}
      $mins{$cond}{$rep} = $min;
      #print "$cond $rep $min\n";
    }
  }

  return \%cnt, \%mins
}
  

=head1 SYNOPSIS

  levelling_lanes.pl    
    
=head1 OPTIONS

=over 4

=item B<-levdir>=I<pathname>

Directory for levelled data.

=item B<-logdir>=<pathname>

Directory for log files.

=item B<-script>=I<pathname>

Path and name of the script C<one_levelling.pl>.

=item B<-queue>=I<string>

Cluster queue name(s), comma delimited.

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back
  