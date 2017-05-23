#!/usr/bin/perl

=head1 NAME

dump2csv.pl - dump bootstrap SQLite data to delimited datafiles

=cut

use strict;
use warnings;

use Getopt::Long qw(:config auto_version);
use Pod::Usage;
use File::Basename;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use DBI;

my $path;
my $out;
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
our $VERSION = '0.1';

GetOptions (
   'path=s'    => \$path,
   'out=s'     => \$out,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid filename.') unless ($path && -d $path);

opendir(my $dh, "$path") or die "ERROR - unable to open directory '$path': $!\nDied";
my @files;
while (my $entry = readdir($dh)) {
   next if ($entry =~ /^\./);  # skip hidden/special files
   push @files, $entry if ($entry =~ /\.db$/);
}
closedir($dh);

printf "Found %d database files\n", scalar @files;

for my $db (@files) {
   my $out = basename($db, ".db");
   $out .= ".csv";
   open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
   print $OUT "BootstrapCount,Gene,WT1,WT2,log2FoldChange,FDR\n";
   my $dbh = DBI->connect("dbi:SQLite:dbname=$path/$db", '', '') or die "ERROR - unable to connect: ", DBI::errstr();
   my $sth = $dbh->prepare("
      SELECT bsid as BScount, f.featureID as gene, WT1, WT2, log2foldchange, significance  
      FROM DEresults de 
      JOIN features f on f.id = de.featureid 
      limit 30"
      ) or die "ERROR - prepare() statement failed: ", $dbh->errstr(); 
   $sth->execute() or die "ERROR - execute() statement failed: ", $dbh->errstr;  
   while (my $row = $sth->fetchrow_arrayref()) {
      print $OUT join(",",@$row),"\n";
   }
   close($OUT);
   $sth->finish();
   $dbh->disconnect();
}

=head1 SYNOPSIS

dump2csv.pl --path <path> [--out <file>] [--verbose|--no-verbose] [--version] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION


=head1 OPTIONS

=over 5

=item B<--in>

Input file.

=item B<--out>

Output filename. [default: STDOUT]

=item B<--version>

Report version information and exit

=item B<--verbose|--no-verbose>

Toggle verbosity. [default:none]

=item B<--debug|--no-debug>

Toggle debugging output. [default:none]

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back

=head1 AUTHOR

Chris Cole <c.cole@dundee.ac.uk>

=head1 COPYRIGHT

Copyright 2017, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut