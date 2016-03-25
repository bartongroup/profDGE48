package BootstrapDB;

=head1 NAME

B<BootstrapDB> - set of routines to create/manipulate Sqlite databases for DE bootstrap

=cut

require Exporter;

use strict;
use warnings;

use DBI;

use PDL;
use PDL::NiceSlice;
use File::Temp qw(tempfile mktemp);

use Time::HiRes qw(gettimeofday tv_interval);

use Tools;

our @ISA = ("Exporter");
our @EXPORT = qw(
  bsConnect
  bsCreateTables
  bsFillFeatures
  bsGetFeatures
  bsAggregate
);

my %tab_name = (
  features => 'features',
  bslogs => 'bslogs',
  bootstraps => 'bootstraps',
  DEresults => 'DEresults'
);

my %tab_sql = (
  features => "CREATE TABLE $tab_name{features} (id INTEGER PRIMARY KEY, featureID TEXT, name TEXT, desc BLOB)",
  bslogs => "CREATE TABLE $tab_name{bslogs} (id INTEGER PRIMARY KEY, log BLOB)",
  bootstraps => "CREATE TABLE $tab_name{bootstraps} (id INTEGER PRIMARY KEY, bsinlogid INTEGER, comments BLOB, deruntime_sec REAL, logid INTEGER, FOREIGN KEY(logid) REFERENCES bslogs(id))",
  DEresults => "CREATE TABLE $tab_name{DEresults} (id INTEGER PRIMARY KEY, featureid INTEGER, bsid INTEGER, log2foldchange REAL, significance REAL, FOREIGN KEY(featureid) REFERENCES features(id), FOREIGN KEY(bsid) REFERENCES bootstraps(id))"
);


###########################################

sub bsConnect
{
  my $dbfile = shift;
  
  my $db = DBI->connect("DBI:SQLite:dbname=$dbfile", '', '', {AutoCommit => 0, RaiseError => 1})  or die   
}

sub bsCreateTables
{
  my $db = shift;
  
  for my $tab (keys %tab_sql)
    {$db->do($tab_sql{$tab}) or die "Failed creating table $tab_name{$tab}\n"}  
}


sub bsFillFeatures
{
  my ($db, $genes) = @_;
  
  my %geneid = ();
  my $insgen = $db->prepare("INSERT INTO $tab_name{features} VALUES (?, ?, ?, ?)") or die "Prepare failed\n";
  for my $geneid (1 .. @$genes)
  {
    my $gene = $genes->[$geneid-1];
    $gene =~ tr/A-Z/a-z/;
    $geneid{$gene} = $geneid;
    $insgen->execute($geneid, $gene, undef, undef) or die "Cannot insert gene $gene\n" 
  }
  $insgen->finish();
  $db->commit();
  
  return %geneid
}

sub bsGetFeatures
{
  my $db = shift;
  
  my $fids = $db->selectall_arrayref("select id, featureID from $tab_name{features}") or die;
  my %featid = ();
  for my $f (@$fids)
  {
    my ($id, $gene) = @$f;
    $featid{$gene} = $id;
  }
  
  return $fids, \%featid;
}


sub bsAggregate
{
  my ($db, $tmpfiles, $commands, $times, $comment, $fccol, $pcol) = @_;
  
  $fccol = 1 unless defined $fccol;
  $pcol = 2 unless defined $pcol;
  $comment ||= '';
  my $nsim = @$tmpfiles;
  
  my ($fids, $geneid) = bsGetFeatures($db);
  
  my $logid = 1;
  $db->do("INSERT INTO $tab_name{bslogs} VALUES ($logid, \"$comment\")") or die "Failed inserting into table $tab_name{bslogs}\n";

  my $insbs = $db->prepare("INSERT INTO $tab_name{bootstraps} VALUES (?, ?, ?, ?, ?)") or die;
  my $insde = $db->prepare("INSERT INTO $tab_name{DEresults} VALUES (?, ?, ?, ?, ?)") or die;

  my $denum = 1;
  for my $bsnum (1 .. $nsim)
  {
    my $tmpfile = $tmpfiles->[$bsnum - 1];
    next unless defined $tmpfile;
    my $cmd = $commands->[$bsnum - 1];
    my $tim = $times->[$bsnum - 1];
    $tim ||= 0;
    
    $cmd = (defined $cmd) ? "Created with $cmd" : '';
    $insbs->execute($bsnum, undef, $cmd, $tim, $logid);
    
    die "Cannot find file $tmpfile\n" unless -e $tmpfile;
    open F, $tmpfile or die;
    while(my $line = <F>)
    {
      chomp $line;
      #my ($gene, $fc, $p, $m1, $m2) = split /\t/, $line;
      my @s = split /\t/, $line;
      my $gene = $s[0];
      my $fc = $s[$fccol];
      my $p = $s[$pcol];
      $gene =~ tr/A-Z/a-z/;
      my $genid = $geneid->{$gene};
      unless(defined $genid)
      {
        #warn "No gene id for $gene\n";
        next;
      }
      $insde->execute($denum, $genid, $bsnum, $fc, $p) or die "Failed inserting into table $tab_name{DEresults}\n";
      $denum++;
    }
    close F;
    #unlink $tmpfile
  }
  $insbs->finish();
  $insde->finish();
  $db->commit();
}


