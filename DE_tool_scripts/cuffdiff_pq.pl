#!/sw/bin/perl

=head1 NAME

Replace p-values (column 3) with q-values (column 4) in cuffdiff output

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use CompBio::Tools;


$| = 1;

my ($infile, $outfile);
my ($help, $man);
GetOptions(
  'infile=s' => \$infile,
  'outfile=s' => \$outfile,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;

Stamp();

die "Need valid -infile\n" unless defined $infile && -e $infile;

open F, $infile or die;
open O, ">$outfile" or die;
while(my $line = <F>)
{
  chomp $line;
  my @s = split /\t/, $line;
  swap($s[2], $s[3]);
  print O join("\t", @s), "\n";
}

sub swap
{
  my $x= $_[0];
  $_[0] = $_[1];
  $_[1] = $x
}


=head1 OPTIONS

=over 5

=item B<-infile>=I<pathname>

Input file.

=item B<-outfile>=I<pathname>

Output file.

=back

=cut
