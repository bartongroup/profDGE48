#!/sw/bin/perl

=head1 NAME

B<make_cuffdiff_db.pl>

=head1 DESCRIPTION

Creates databases from cuffdiff bootstrap results. You need to run C<grid_launcher_powertest_cuffdiff.pl> first. 

=cut


use strict;
use warnings;

use PDL;
use PDL::NiceSlice;

use Cwd;
use Getopt::Long;
use Pod::Usage;

use GRNASeq;
use CompBio::Tools;
use BootstrapDB;

$| = 1;
Stamp();

my $cwd = getcwd; 
my $cuffdir = "$cwd/tmp_cuffdir";
my $dbdir = "$cwd/de_tests_db";
my $genlist;
my $sig = 'q';
my ($help, $man);
GetOptions(
  'cuffdir=s' => \$cuffdir,
  'dbdir=s' => \$dbdir,
  'genlist=s' => \$genlist,
  'sig=s' => \$sig,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;

die "Need valid -genlist file\n" unless defined $genlist && -e $genlist;
my @genes = ReadGeneList($genlist);
opendir D, $cuffdir or die "Cannot open directory $cuffdir\n";
my @f = readdir D;

die "-sig must be p or q\n" unless $sig eq 'p' || $sig eq 'q';
my $pcol = ($sig eq 'p') ? 2 : 3;

print "Scanning cuffdiff output...";
my @tfiles = ();
my @cmd = ();
my $minrep = 100;
my $maxrep = -1;
for my $d (@f)
{
  my $dir = "$cuffdir/$d";
  next unless -d $dir;
  if($d =~ /cuffdiff_out_r(\d+)_i(\d+)/)
  {
     my $nrep = $1;
     my $iter = $2;
     $d =~ /cuffdiff_out_(.*)/;
     my $num = $1;
     my $tfile = "$cuffdir/cuffdiff_$num.tsv";
     if(-e $tfile)
     {
       my $cfile = "$dir/README";
       open F, $cfile or die "Cannot open $cfile\n";
       my $cmd = <F>;
       chomp $cmd;
       close F; 
       
       $tfiles[$nrep][$iter-1] = $tfile;
       $cmd[$nrep][$iter-1] = $cmd;
       
       $minrep = $nrep if $nrep < $minrep;
       $maxrep = $nrep if $nrep > $maxrep;
     }
     #else
     #{print "No $tfile\n"}
     
     #print "$nrep $iter $tfile\n";
  }
}
print " done\n\n";

for my $nrep (1 .. @tfiles-1)
{
  my $sum = 0;
  if(defined $tfiles[$nrep])
  {
    my @s = @{$tfiles[$nrep]};
    for my $s (@s)
      {$sum++ if defined $s}
  }
  print "$nrep: $sum\n";
}
print "\n";

print "Building db files\n\n";
for my $nrep ($minrep .. $maxrep)
{
  next unless defined $tfiles[$nrep];
  print "nrep = $nrep\n";

  my @t = @{$tfiles[$nrep]};
  my @c = @{$cmd[$nrep]};

  my $q = ($sig eq 'q') ? 'q' : '';
  my $dbfile = "$dbdir/de_cuffdiff${q}_rep$nrep.db";
  print "  Creating database $dbfile...";
  unlink $dbfile;
  my $db = bsConnect($dbfile);
  bsCreateTables($db);
  bsFillFeatures($db, \@genes);
  print " done\n";

  print "  Aggregating results...";
  bsAggregate($db, \@t, \@c, undef, undef, undef, $pcol);
  print " done\n";
}



=head1 SYNOPSIS

  make_cuffdiff_db.pl -cuffdir=cuffdiff_powertests -dbdir=de_tests_db -genlist=genes.lst   

=head1 OPTIONS

=over 5

=item B<-cuffdir>=I<path>

Path to the directory where cuffdiff bootstrap results are stored.

=item B<-dbdir>=I<path>

Path to the directory where all results are saved. These are Sqlite databases, one file per number of replicates.

=item B<-genlist>=I<file>

File with gene names in the first column.

=item B<-sig>=[p|q]

Whether to use p-value or q-value from the cuffdiff's output. Default is C<p>.

=item B<-help>

Brief help.

=item B<-man>

Full manpage of program.

=back

=cut