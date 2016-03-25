#!/sw/bin/perl

=head1 NAME

B<combine replicates>

=head1 DESCRIPTION

Combine counts from individual condition/replicate/lane files into multicolumn files. Genes in all files B<must> be in the same order.

Warning: replicate and lane number and file name templates hardcoded in the script (sorry...)

Usage:

  combine_replicates.pl -indir=<input dir> -outdir=<output dir> -bioreps
  
-bioreps tells the script to do biological replicates instead of lanes.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;

use CompBio::Tools;
use GRNASeq;

$| = 1;
Stamp();

my @conds = qw(WT Snf2);
my @reps = (1 .. 48);
my @lanes = (1..7);
my $post = 'tophat2.0.5_genes';

my $filter = '(no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique)';

my ($help, $man);
my ($indir, $outdir);
my ($bioreps, $fpkm);
GetOptions(
  'indir=s' => \$indir,
  'outdir=s' => \$outdir,
  'post=s' => \$post,
  bioreps => \$bioreps,
  fpkm => \$fpkm,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;

die "Need -indir\n" unless defined $indir;
die "Need -outdir\n" unless defined $outdir;

#$indir = "$topdir/analysis/$indir";
#$outdir = "$topdir/analysis/$outdir";

die "Cannot find $indir\n" unless -d $indir;
die "Cannot find $outdir\n" unless -d $outdir;

my $N = scalar @conds * scalar @reps * scalar @lanes;

my $ext = ($fpkm) ? 'fpkm' : 'tsv';
my $form_replane = "%s/%s_rep%02d_lane%1d%s.$ext";
my $form_rep = "%s/%s_rep%02d_allLanes_%s.$ext";

my ($col, $skip);
if($fpkm)
{
  $col = 9;            # column with counts
  $skip = 1;           # lines to skip at the start of the file
}
else
{
  $col = 1;
  $skip = 0;
}


if($bioreps)
{
  CombineBioReplicates();
}
else
{
  CombineLanes();
  CombineLaneReplicates();
}

#####################################

sub CombineLanes
{
  my $firstfile = sprintf $form_replane, $indir, $conds[0], 1, 1, $post;
  die "Cannot find $firstfile\n" unless -e $firstfile;
  my @genes = ReadGeneList($firstfile, $skip);
  
  my $cnt = 1 ;
  print "Combining lanes...       ";
  for my $cond (@conds)
  {
    for my $rep (@reps)
    {
      my @d = ();
      for my $lane (@lanes)
      {
        Percent1(($cnt++)/$N);
        my $file = sprintf $form_replane, $indir, $cond, $rep, $lane, $post;
        my $x = QuickReadCountFile($file);
        push @d, $x 
      }
      my $d = pdl(@d)->transpose();
      my ($cols, $rows) = dims $d;
      my $out = sprintf "%s/%s_rep%02d_raw.tsv", $outdir, $cond, $rep;
      my $no = sprintf "%s/%s_rep%02d_nogene.tsv", $outdir, $cond, $rep;
      local (*F, *N);
      open F, ">$out" or die;
      open N, ">$no" or die;
      for my $i (0 .. $rows - 1)
      {
        my $row = $genes[$i] . "\t" . join("\t", list($d(,$i;-)));
        if($genes[$i] =~  /$filter/)
          {print N "$row\n"}
        else
          {print F "$row\n"}
      }
      close F;
    }
  }
  print "\b\b\b\b\b\b\b done   \n";
}

#####################################

sub CombineLaneReplicates
{
  my $firstfile = sprintf $form_replane, $indir, $conds[0], 1, 1, $post;
  die "Cannot find $firstfile\n" unless -e $firstfile;
  my @genes = ReadGeneList($firstfile, $skip);
  
  my $cnt = 1;
  print "Combining replicates...       ";
  for my $cond (@conds)
  {
    my @c = ();
    for my $rep (@reps)
    {
      my @d = ();
      for my $lane (@lanes)
      {
        Percent1(($cnt++)/$N);
        my $file = sprintf $form_replane, $indir, $cond, $rep, $lane, $post;
        my $x = QuickReadCountFile($file);
        push @d, $x 
      }
      my $d = pdl(@d)->transpose();
      my $s = sumover($d);
      push @c, $s;
    }
    my $c = pdl(@c)->transpose();
    my ($cols, $rows) = dims $c;
    my $out = sprintf "%s/%s_raw.tsv", $outdir, $cond;
    my $no = sprintf "%s/%s_nogene.tsv", $outdir, $cond;
    local (*F, *N);
    open F, ">$out" or die;
    open N, ">$no" or die;
    open F, ">$out" or die;
    for my $i (0 .. $rows - 1)
    {
      my $row = $genes[$i] . "\t" . join("\t", list($c(,$i;-)));
      if($genes[$i] =~  /$filter/)
        {print N "$row\n"}
      else
        {print F "$row\n"}
    }
    close F;
  }
  print "\b\b\b\b\b\b\b done   \n";
}

#####################################

sub CombineBioReplicates
{
  my $firstfile = sprintf $form_rep, $indir, $conds[0], 1, $post;
  die "Cannot find $firstfile\n" unless -e $firstfile;
  my @genes = ReadGeneList($firstfile, $skip);

  my $b = "\b" x 7;
  my $cnt = 1 ;
  print "Combining biological replicates...       ";
  for my $cond (@conds)
  {
    my @c = ();
    for my $rep (@reps)
    {
      #Percent1(($cnt++)/$N);
      printf "$b%4s %2d", $cond, $rep;
      my $file = sprintf $form_rep, $indir, $cond, $rep, $post;
      my $x = QuickReadCountFile($file);
      push @c, $x 
    }
    my $c = pdl(@c)->transpose();
    my ($cols, $rows) = dims $c;
    my $out = sprintf "%s/%s_raw.tsv", $outdir, $cond;
    my $no = sprintf "%s/%s_nogene.tsv", $outdir, $cond;
    local (*F, *N);
    open F, ">$out" or die;
    open N, ">$no" or die;
    for my $i (0 .. $rows - 1)
    {
      my $row = $genes[$i] . "\t" . join("\t", list($c(,$i;-)));
      if($genes[$i] =~  /$filter/)
        {print N "$row\n"}
      else
        {print F "$row\n"}
    }
    close F;
  }
  print "${b} done   \n";
}

#####################################
#
# Careful: we assume genes in all files are in the same order!
#
#
sub QuickReadCountFile
{
  my ($file) = @_;

  die "Missing file $file\n" unless -e $file;
  my $x = rcols $file, $col, {LINES=>"$skip:-1"}
}

