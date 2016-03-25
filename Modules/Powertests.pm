package Powertests;

=head1 NAME

B<Powertests>

=cut

require Exporter;

use strict;
use warnings;

use PDL;
use PDL::NiceSlice;

use Tools;

our @ISA = ("Exporter");
our @EXPORT = qw(
  rerr
  GetTotData
  GetFCData
  GetFCOneData
  CollectFC
);


##########################################

=item B<rerr>

  my ($r, $sr) = rerr($x, $y, $sx, $sy);

Calculates $x / ($x + $y) and its error.

=cut


sub rerr
{
  my ($x, $y, $sx, $sy) = @_;
  
  my $S = $x + $y;
  my $r = $x / $S;
  my $sr = sqrt($sx**2 * (1 / $S - $x / $S**2)**2 + ($sy * $x)**2 / $S**4);
  
  #$r->where($S == 0) .= 0;
  #$sr->where($S == 0) .= 0;
  
  return $r, $sr;
}

##########################################

=item B<GetTotData>

  my ($nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetTotData($file)

Reads data from powertest file.

=cut

sub GetTotData
{
  my $file = shift;
  
  my ($nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = rcols $file, {LINES=>'1:-1'}; 
#  my ($fpr, $sfpr) = rerr($fp, $tn, $sfp, $stn);
#  my ($fnr, $sfnr) = rerr($fn, $tp, $sfn, $stp);
  
  return $nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn
}

##########################################

=item B<GetFCData>

  my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetTotData($file)

Reads data from powertest FC file.

=cut

sub GetFCData
{
  my $file = shift;
  
  my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = rcols $file, {LINES=>'1:-1'};
  
  return $nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn
}


##########################################

=item B<GetFCOneData>

  my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetTotData($file)

Reads data from powertest FC file and combines positive and negative FC.

=cut

sub GetFCOneData
#
# Read FC data and sum both tails
# The result is rates for abs(FC) > fclim
#
{
  my $file = shift;
  
  my ($nrep_, $fclim_, $tp_, $tn_, $fp_, $fn_, $stp_, $stn_, $sfp_, $sfn_) = rcols $file, {LINES=>'1:-1'};
  
  my @nrep = ();
  my @fclim = ();
  my @tp = ();
  my @tn = ();
  my @fp = ();
  my @fn = ();
  my @stp = ();
  my @stn = ();
  my @sfp = ();
  my @sfn = ();
  
  my @F = CollectFC($fclim_);
  for my $n (2 .. max($nrep_))
  {
    for my $f (@F)
    {
      next if $f < 0;
      my $sel = which((abs($fclim_) == $f) & ($nrep_ == $n));
      next unless nelem($sel) == 2;
      push @nrep, $n;
      push @fclim, $f;
      push @tp, sum($tp_($sel));
      push @fp, sum($fp_($sel));
      push @tn, sum($tn_($sel));
      push @fn, sum($fn_($sel));
      push @stp, sqrt(sum($stp_($sel)**2));
      push @sfp, sqrt(sum($sfp_($sel)**2));
      push @stn, sqrt(sum($stn_($sel)**2));
      push @sfn, sqrt(sum($sfn_($sel)**2));
    }
  }
  return pdl(@nrep), pdl(@fclim), pdl(@tp), pdl(@tn), pdl(@fp), pdl(@fn), pdl(@stp), pdl(@stn), pdl(@sfp), pdl(@sfn)
}

##########################################

=item B<CollectFC>

  my @F = CollectFC($fclim);

Collects unique fold changes from a piddle.

=cut



sub CollectFC
{
  my $fclim = shift;
  
  my $sfclim = qsort $fclim;
  my @F = ();
  my $prev = -1e16;
  for my $i (0 .. nelem($sfclim) - 1)
  {
    my $f = at($sfclim, $i);
    if($f != $prev)
    {
      push @F, $f;
      $prev = $f
    }
  }
  return @F
}



1;

