#!/sw/bin/perl -w

=head1 NAME

B<one_levelling.pl> - script used by C<levelling_lanes.pl>. Not to be run on its own.  

=cut

#$ -e /cluster/gjb_lab/mgierlinski/projects/grna/NOBACK/gridsim -o /cluster/gjb_lab/mgierlinski/projects/grna/NOBACK/gridsim

use strict;
use warnings;

use PDL;
use PDL::NiceSlice;
use File::Temp qw(tempfile mktemp);

use CompBio::Tools;

use Getopt::Long;
use Pod::Usage;

my ($infile, $outfile, $count, $min);
my ($help, $man);
GetOptions (
  'infile=s' => \$infile,
  'outfile=s' => \$outfile,
  'count=i' => \$count,
  'min=i' => \$min,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;



print "Running levelling\n";
print "Infile = $infile\n";
print "Outfile = $outfile\n";
print "Count = $count\n";
print "Min = $min\n";


Level($infile, $outfile, $count, $min);


####################################################

sub Level
{
  my ($infile, $outfile, $nreads, $ntarget) = @_;
  
  print "randomizing... ";
  my $r = random($nreads);    # random vector
  my $perm = qsorti $r;       # random permutation of $nreads indices
  my $sel = qsort $perm(0:$ntarget-1);   # $ntarget random indices
  
  my $tmpfile = mktemp("./levelXXXXXX");
    
  print "processing...       ";
  local (*IN, *OUT);
  open IN, "gunzip -c $infile |" or die;
  open OUT, ">$tmpfile" or die;
  my $i = 0;
  my $n = 0;
  while($i < $ntarget - 1)
  {
    Percent1($i / ($ntarget-1)) if $i % 10000 == 0;
    my $read = GetRead(*IN);
    if(at($sel, $i) == $n)
    {
      #print "$i $n\n";
      print OUT "$read";
      $i++
    }
    $n++
  }
  
  print "\b\b\b\b\b\b\b gzipping...";
  `gzip -c $tmpfile > $outfile`;
  unlink $tmpfile; 
  
  print "\n";
  if(-e $outfile && !-e $tmpfile)
    {print "Probably all done\n"}
  else
    {die "Something went wrong\n"}
}

sub GetRead
{
  my $z = shift;
  my $r = <$z>;
  die "Unexpected end of file\n" unless defined $r;
  for my $i (2 .. 4)
    {$r .= <$z>}
  return $r;
}


=head1 OPTIONS

=over 4

=item B<-infile>=I<pathname>

Input file.

=item B<-outfile>=<pathname>

Output file.

=item B<-count>=I<number>

Number of read counts.

=item B<-min>=I<number>

Target number of reads.



=back
  
