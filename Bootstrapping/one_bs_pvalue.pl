#!/sw/bin/perl

=head1 NAME

B<one_bs_pvalue.pl>

=head1 DESCRIPTION

To be used by B<combine_bs_pvalues.pl>. This script is not supposed to be run directly.

=cut


use strict;
use warnings;

use DBI;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;

use CompBio::Tools;
use CompBio::Statistics;

use Getopt::Long;

$| = 1;

my ($dbfile, $outfile, $logit);
my ($help, $man);
GetOptions(
  'dbfile=s' => \$dbfile,
  'outfile=s' => \$outfile,
  logit => \$logit,
  help => \$help,
  man => \$man, 
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;


die "Need dbfile\n" unless defined $dbfile && -e $dbfile;
die "Need outfile\n" unless defined $outfile;

print "Combining bootstrap results for $dbfile\n";
CreateMedianPvalues($dbfile, $outfile);

#################################

sub CreateMedianPvalues
{
  my ($dbfile, $outfile) = @_;
  
  #my @columns = qw(featureid bsid log2foldchange significance meanc1 meanc2);
  #my $columns = join ",", @columns;
  
  my $db = DBI->connect("DBI:SQLite:dbname=$dbfile", '', '', {AutoCommit => 0, RaiseError => 1})  or die "\nCouldn't open local database: " . DBI->errstr;
  
  # get geneid, gene name
  my $fids = $db->selectall_arrayref("select id, featureID from features") or die;
  
  # DEresults structure, find additional columns
  #my $tcols = $db->selectall_arrayref("pragma table_info(DEresults)");
  
  
  # get geneid, pvalue
  my $rows = $db->selectall_arrayref("select featureid, log2foldchange, significance from DEresults") or die;
  my $dat = pdl(@$rows);
  my $featureid = $dat(0,;-);
  my $fc = $dat(1,;-);
  if($logit)
  {
    $fc->where($fc == 0) .= 1;  # SamSeq reports fc=0 for some weak genes
    $fc = log($fc)/log(2);
  }
  my $pvalue = $dat(2,;-);

  local *F;
  open F, ">$outfile" or die "\nCannot open output file $outfile\n";

  # cycle through all feature ids (genes) and find geometric means
  my $N = @$fids;
  my $cnt = 0;
  for my $r (@$fids)
  {
    Percent1(($cnt++)/$N);
    my ($id, $gene) = @$r;
    my $pid = $pvalue->where($featureid == $id);
    my $fcid = $fc->where($featureid == $id);
    next if nelem($pid) == 0;
    my $meanp = exp(mean(log($pid)));
    my $meanfc = mean($fcid);
    my $medp = median($pid);
    #print F "$gene\t$meanfc\t$meanp\n";
    print F "$gene\t$meanfc\t$medp\t$meanp\n";  # median and geom mean
  }
  print "\n";
}

=head1 OPTIONS

=over 5

=item B<-dbfile>=I<pathname>

Database file to process.

=item B<-outfile>=I<pathname>

Output file.


=item B<-logit>

If this option is specified, log2  will be calculated of all fold-changes. However, we expect .db files to contain log2-fold-changes, so this is probably not necessary.


=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut
