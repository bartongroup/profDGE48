#!/sw/bin/perl

=head1 NAME

B<one_bs_powerstats.pl>

=head1 DESCRIPTION

To be used by C<make_powerstats_db.pl>. Not to be run on its own.

=cut


use strict;
use warnings;

use DBI;

use PDL;
use PDL::NiceSlice;

use Tools;
use Stats;
use GRNASeq;

use Getopt::Long;
use Pod::Usage;

$| = 1;

my $version = 1.1;

my $created = "Created with $0 " . join(" ", @ARGV) . "\n";

my $alpha = 0.05;
my $multicor = 'bh';
my ($truecor, $dbcor);
my ($dbfile, $truefile, $totfile, $fcfile1, $fcfile2, $rocfile, $fcsig, $reffcfile, $nsigfile, $sigpropfile);
my $genlist = 'genlist.tsv';
my $colfc = 2;
my $colp = 3;
my ($help, $man);

GetOptions(
  'genlist=s' => \$genlist,
  'dbfile=s' => \$dbfile,
  'totfile=s' => \$totfile,
  'fcfile1=s' => \$fcfile1,
  'fcfile2=s' => \$fcfile2,
  'rocfile=s' => \$rocfile,
  'truefile=s' => \$truefile,
  'nsigfile=s' => \$nsigfile,
  'sigpropfile=s' => \$sigpropfile,
  'multicor=s' => \$multicor,
  'truecor=s' => \$truecor,
  'dbcor=s' => \$dbcor,
  'alpha=f' => \$alpha,
  'fcsig=f' => \$fcsig,
  'reffcfile=s' => \$reffcfile,
  'colp=i' => \$colp,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;


die "Need dbfile\n" unless defined $dbfile && -e $dbfile;
die "Need totfile\n" unless defined $totfile;
die "Need nsigfile\n" unless defined $nsigfile;
die "Need sigpropfile\n" unless defined $sigpropfile;
die "Need fcfile1\n" unless defined $fcfile1;
die "Need fcfile2\n" unless defined $fcfile2;
die "Need rocfile\n" unless defined $rocfile;
die "Need truefile\n" unless defined $truefile;
die "Truefile $truefile does not exist\n" unless  -e $truefile;

print $created;
print "Calculating power statistics for $dbfile using true test from $truefile\n";


my $db = DBI->connect("DBI:SQLite:dbname=$dbfile", '', '', {AutoCommit => 0, RaiseError => 1})  or die "\nCouldn't open local database: " . DBI->errstr;
  
my @genes = ReadGeneList($genlist);
  
my ($f2i, $ptr, $fctr, $sel, $selgenes) = TrueFileIndex($db, $genlist, $truefile);
my $N = nelem($ptr);
# f2i is a cross-reference database's feature id -> gene id from genlist

my $reffc;
if($reffcfile)
{
  my ($reff2i, $pref, $fref) = TrueFileIndex($db, $genlist, $reffcfile);
  $reffc = $fref;
}

# warning: here we assume that "true" fold change does not need sign correction!

my $ptrsel = $ptr($sel);
my $fctrsel = $fctr($sel);
$fctrsel = $reffc($sel) if defined $reffcfile;

#my $test = sum(($fctr($sel) - $reffc($sel))**2);
#print "TEST=$test\n";

#wcols $fctr($sel), $reffc($sel);die;


my $rows = $db->selectall_arrayref("select featureid, bsid, log2foldchange, significance from DEresults where significance!='NA'") or die;
my $dat = pdl(@$rows);
my $featureid = $dat(0,;-);
my $bsid = $dat(1,;-);
my $fc = $dat(2,;-);
my $pvalue = $dat(3,;-);

$fc *= $fcsig;

my $nbs = max($bsid);
my @fcstat1 = ();
my @fcstat2 = ();
my @totstat = ();
my @rocstat = ();
my @nsig = ();
my @sigprop = ();  # proportion of runs showing significance (per gene)


my ($fcmin, $fcmax, $fcstep) = (-4, 4, 0.05);
my $fcn1 = ($fcmax - $fcmin) / $fcstep + 1;
my $fcn2 = ($fcmax - 0) / $fcstep + 1;

my @fcsel = (0, 0.5, 1, 2);
my $fcnsel = @fcsel;

#
# Bootstrap loop
#
#$nbs = 7;
for my $bs (1 .. $nbs)
{
  Percent1($bs/$nbs);
  my $bs_fid = $featureid->where($bsid == $bs);
  my $bs_p = $pvalue->where($bsid == $bs);
  my $bs_fc = $fc->where($bsid == $bs);

  my $p = zeroes($N) + 1;
  my $fc = zeroes($N);

  $p($f2i($bs_fid)) .= $bs_p;   # reindexing to gene lists's gene order
  $fc($f2i($bs_fid)) .= $bs_fc;
  
  # select only genes present in truefile
  my $psel = $p($sel);
  my $fcsel = $fc($sel);
  $fcsel = $reffc($sel) if defined $reffcfile;
  
  # now ptrsel and psel have the same order and selection
  # so we can compare them directly
  
  #print "N = ", nelem($bs_p), " ",nelem($p), " ", nelem($ptr), "\n";
  #wcols $bs_fid(0:12), $bs_p(0:12);
  #print "f2i=", $f2i(0:12), "\n";
  #my $bb = $f2i($bs_fid);
  #print "f2i(bs_fid)=", $bb(0:12), "\n";
  #wcols $sel(0:20), $psel(0:20), $ptrsel(0:20);
  #die;

  # significant proportion
  my $lim = MulticorLimit($psel, $dbcor, $multicor, $alpha);
  my $sig = ($dbcor =~ /llim/) ? ($psel >= $lim) : ($psel <= $lim);
  push @sigprop, $sig;

  # power test
  my ($TP, $TN, $FP, $FN, $s1, $s2) = PowerTest($ptrsel, $psel, $fctrsel, $fcsel, $multicor, $truecor, $dbcor, 0, $alpha);
  $totstat[$bs - 1] = [$TP, $TN, $FP, $FN];
  #$nsig[$bs - 1] = $TP + $FP;  # all significant genes
  $nsig[$bs - 1] = ($TP + $FP) / ($TP + $TN + $FP + $FN);  # fraction significant genes
  
#  if($bs == 1)
#  {
#    for my $i (0 .. nelem($sel) - 1)
#    {
#      my $gene = $selgenes->[$i];
#      my $ss1 = at($s1, $i);
#      my $ss2 = at($s2, $i);
#      print "$gene    $ss1  $ss2\n" if $ss2 > $ss1;
#    }
#    die;
#  }
  
  #my $f = 0;
  #for (my $fclim = $fcmin; $fclim <= $fcmax; $fclim += $fcstep)
  
  # power tests for a range of fold-change limits
  # one-sided tests are for fc < fclim, two-sided tests are for fc >= fclim
  
  # one-sided test:
  for my $f (0 .. $fcn1 - 1)
  {
    my $fclim = $fcmin + $fcstep * $f;
    my ($TP, $TN, $FP, $FN) = PowerTest($ptrsel, $psel, $fctrsel, $fcsel, $multicor, $truecor, $dbcor, $fclim, $alpha, 1);
    $fcstat1[$bs - 1][$f] = [$TP, $TN, $FP, $FN];
  }
  
  # two-sided test
  for my $f (0 .. $fcn2 - 1)
  {
    my $fclim = 0 + $fcstep * $f;
    my ($TP, $TN, $FP, $FN) = PowerTest($ptrsel, $psel, $fctrsel, $fcsel, $multicor, $truecor, $dbcor, $fclim, $alpha);
    $fcstat2[$bs - 1][$f] = [$TP, $TN, $FP, $FN];
  }
  
  # ROC curves
  for my $f (0 .. $fcnsel - 1)
  {
    my $limit = 1e-5;
    my @fpr = ();
    my @tpr = ();
    #my @tnr = ();
    while($limit < 1)
    {  
      #my $fclim = 0 + $fcstep * $f;
      my $fclim = $fcsel[$f];
      my ($TP, $TN, $FP, $FN) = PowerTest($ptrsel, $psel, $fctrsel, $fcsel, $multicor, $truecor, $dbcor, $fclim, $limit);
      my $FP_rate = ($FP + $TN) > 0 ? $FP / ($FP + $TN) : 1;
      my $TP_rate = ($TP + $FN) > 0 ? $TP / ($TP + $FN) : 0;
      #my $TN_rate = ($TN + $FP) > 0 ? $TN / ($TN + $FP) : 0;
      last if $FP + $TN == 0;
      #push @tnr, $TN_rate;
      push @tpr, $TP_rate;
      push @fpr, $FP_rate;
      $limit *= 1.1;
    }
    $rocstat[$bs - 1][$f] = [pdl(@fpr), pdl(@tpr)];
  }
}
print "\n";

my $eform = "%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n";


# processing num significant genes

open SIG, ">$nsigfile" or die;
print SIG join("\t", @nsig), "\n";
close SIG;

#processing significant proportion

my $sig = pdl(@sigprop)->transpose();
my $sigprop = sumover($sig) / $nbs;
open PROP, ">$sigpropfile" or die;
for my $i (0 .. nelem($sigprop) - 1)
{
  my $gene = $genes[at($sel, $i)];
  my $fc = at($fctrsel, $i);
  my $prop = at($sigprop, $i);
  print PROP "$gene\t$fc\t$prop\n";
}
close PROP;

# processing fc bootstraps
# $fcstat has 3 dimensions: 4 x $fcn x $nbs
# calculating statistics along bootstraps 

open FC1, ">$fcfile1" or die;
my $fcstat1 = pdl(@fcstat1)->xchg(0,2); # select bootstraps as 1st dimension
my ($mfc1, $sfc1) = statsover($fcstat1);  # stats over bootstraps
for my $f (0 .. $fcn1 - 1)
{
  my $fclim = $fcmin + $fcstep * $f;
  printf FC1 "%4.2f\t", $fclim;
  printf FC1 $eform, list($mfc1($f,;-)), list($sfc1($f,;-));
}
close FC1;

open FC2, ">$fcfile2" or die;
my $fcstat2 = pdl(@fcstat2)->xchg(0,2); # select bootstraps as 1st dimension
my ($mfc2, $sfc2) = statsover($fcstat2);  # stats over bootstraps
for my $f (0 .. $fcn2 - 1)
{
  my $fclim = 0 + $fcstep * $f;
  printf FC2 "%4.2f\t", $fclim;
  printf FC2 $eform, list($mfc2($f,;-)), list($sfc2($f,;-));
}
close FC2;


# processing total bootstraps
# $totstat has 2 dimensions: 4 x $nbs
# calculating statistics along bootstraps

open TOT, ">$totfile" or die;
my $totstat = pdl(@totstat)->xchg(0,1);  # select bootstraps as 1st dimension
my ($mtot, $stot) = statsover($totstat);
printf TOT $eform, list($mtot), list($stot);
close TOT;

# processing ROC bootstraps
# $rocstats has 4 dimensions: alpha_range x 2 x $fcn x $nbs

my $rocstat = pdl(@rocstat)->xchg(0,3);
my ($mroc, $sroc) = statsover($rocstat);
my ($n1, $n2, $n3) = dims($mroc);
open ROC, ">$rocfile" or die;
for my $f (0 .. $fcnsel - 1)
{
  #my $fclim = 0 + $fcstep * $f;
  my $fclim = $fcsel[$f];
  for my $i (0 .. $n3 - 1)
  {
    printf ROC "%4.2f\t%.3g\t%.3g\t%.3g\t%.3g\n", $fclim, at($mroc, 0, $f, $i), at($mroc, 1,$f, $i), at($sroc, 0, $f, $i), at($sroc, 1,$f, $i)
  }
}
print "Created $totfile, $fcfile1, $fcfile2, $rocfile, $nsigfile and $sigpropfile\n";

#################################################

sub TrueFileIndex
{
  my ($db, $genlist, $file) = @_;
  
  # build gene->index hash for all possible genes
  my @allgenes = ReadGeneList($genlist);
  my $N = @allgenes;
  my %g2i = ();
  for my $i (0 .. $N - 1)
    {$g2i{$allgenes[$i]} = $i}
  
  my $fids = $db->selectall_arrayref("select id, featureID from features") or die;
  print scalar @$fids, " rows read\n";
  my @f2i = ();
  for my $f (@$fids)
  {
    my ($fid, $gene) = @$f;
    $gene =~ s/\"//g;
    $gene =~ tr/A-Z/a-z/;
    if(defined $g2i{$gene})
      {$f2i[$fid] = $g2i{$gene}}
    else
      {die "Unidentified gene $gene in features in $dbfile\n"}
  }
      
  local *F;
  my $ptr = zeroes($N) + 1;
  my $fctr = zeroes($N);
  my $is = zeroes($N);
  my @selgenes = ();
  
  open F, $file or die "Cannot open file $file\n";
  my $i=0;
  while(my $line = <F>)
  {
    chomp $line;
    next if $line =~ /^#/;
    next if $line =~ /^gene/i;
    
    my @s = split /\t/, $line;;
    my $gene = $s[0];
    my $fc = $s[$colfc - 1];
    my $p = $s[$colp - 1];
    
    #my ($gene, $fc, $p) = split /\t/, $line;
    
    die "$file: $gene undefined p in line $i $line and file $file\n" unless defined $p;
    $i++;
    
    $p = 1 if $p =~ /NA/;
    $fc = 0 if $fc =~ /NA/;
    $gene =~ s/\"//g;
    $gene =~ tr/A-Z/a-z/;

    my $idx = $g2i{$gene};
    if(defined $idx)
    {
      $ptr($idx) .= $p;
      $fctr($idx) .= $fc; 
      $is($idx) .= 1;
      push @selgenes, $gene;
    }
    else
      {die "Unidentified gene $gene in file $truefile\n"}
  }
  my $sel = which($is);
  
  return pdl(@f2i), $ptr, $fctr, $sel, \@selgenes
}

=head1 OPTIONS

=over 5

=item B<-genlist>=I<pathname>

File with list of genes.

=item B<-dbfile>=I<pathname>

Input sqlite .db file.

=item B<-totfile, -fcfile1, -fcfile2, -rocfile, -truefile, -nsigfile, -sigpropfile>

Various temporary files to communicate with the parent script.

=item B<-multicor>=I<string>

Multiple test correction.

=item B<-truecor>=I<string>

Multiple test correction used in the 'true' file.

=item B<-alpha>=I<number>

Significance limit. Default value is 0.05.

=item B<-fcsig>=I<-1|1>

Sign of log-fold-change.

=item B<-reffcfile>=I<pathname>

The same as in the parent script.

=item B<-colp>=I<number>

Column with p-values.

=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut
  