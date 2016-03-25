#!/sw/bin/perl

=head1 NAME

B<get_gene_info.pl>

=head1 DESCRIPTION

Get gene descriptions from Ensembl and store them in a tab-separated output file with three columns: gene id, gene name and description. The default name of the output file is C<gene_descriptions.tsv> in the current directory.

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use CompBio::Tools;

#use lib '/opt/perl/bioperl-live';
use lib '/sw/opt/ensembl-api/68/ensembl/modules';
use Bio::EnsEMBL::Registry;

$|=1;
Stamp();

my $genlist;
my ($help, $man);
GetOptions(
  'genlist=s' => \$genlist,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help || !defined $genlist;

my $description_file = 'gene_descriptions.tsv';

Bio::EnsEMBL::Registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
) or die;
my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "Saccharomyces cerevisiae", "core", "gene" );

open IN, $genlist or die;
open OUT, ">$description_file" or die;
while(my $line = <IN>)
{
  chomp $line;
  next if $line =~ /^#/;
  my ($gene) = split " ", $line;
  
   my $g = $gene_adaptor->fetch_by_stable_id($gene);
   my $desc = $g->description();
   my $name = $g->external_name();
   $name ||= $gene;
   $desc ||= 'N/A';
   print OUT "$gene\t$name\t$desc\n";
}

=head1 SYNOPSIS

  get_gene_info -genlist=WT_raw.tsv

head1 OPTIONS

=over 4

=item B<-genlist>

Gene list, where first column contains gene names. Other columns are ignored.

=back
