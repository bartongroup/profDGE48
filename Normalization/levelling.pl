#!/sw/bin/perl -w

=head1 NAME

B<levelling.pl>  

Equal count (levelling) normalization of fastq files. No options, just run it. Read count files have to be prepared in advance using C<count_fastq_reads.pl>.

=cut

use strict;
use warnings;

use PDL;
use PDL::NiceSlice;

use CompBio::Tools;

use GRNASeq;

$| = 1;

our $topdir;
my $fastqdir = "$topdir/raw";
my @conds = qw(WT Snf2);

my ($help, $man);
GetOptions(
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;



my ($cnt, $min) = ReadCountFiles();
for my $cond (keys %$cnt)
{
  for my $rep (sort {$a <=> $b} keys %{$cnt->{$cond}})
  {
    my ($infile, $count) = @{$cnt->{$cond}{$rep}};
    my $outfile;
    if($infile =~ /(.*)\.fastq\.gz/)
      {$outfile = "${1}_levelled.fastq.gz"}
    else
      {die "Unrecognized file name $infile\n"}
    print "$cond $rep $count  ";
    Level($infile, $outfile, $count, $min);
    print "\n";
  }
}


####################################################

sub Level
{
  my ($infile, $outfile, $nreads, $ntarget) = @_;
  
  print "randomizing... ";
  my $r = random($nreads);    # random vector
  my $perm = qsorti $r;       # random permutation of $nreads indices
  my $sel = qsort $perm(0:$ntarget-1);   # $ntarget random indices
  
  my $tmpout = "./leveout$$.fastq";
    
  print "processing...       ";
  local (*IN, *OUT);
  open IN, "gunzip -c $infile |" or die;
  open OUT, ">$tmpout" or die;
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
  `gzip -c $tmpout > $outfile`;
  unlink $tmpout; 
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


sub ReadCountFiles
{
  
  my %cnt = ();
  my $min = 1e32;
  for my $cond (@conds)
  {
    local *F;
    my $cfile = "${cond}_readcount.dat";
    open F, $cfile or die;
    while(my $line = <F>)
    {
      chomp $line;
      my ($file, $rep, $count) = split " ", $line;
      $cnt{$cond}{$rep} = [$file, $count];
      $min = $count if $count < $min
    }
  }
  return \%cnt, $min
}
  
