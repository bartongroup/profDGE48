#!/sw/bin/perl

=head1 NAME

run_cuffdiff.pl - wrapper around cuffdiff 

=cut

use strict;
use warnings;

$|=1;

use Getopt::Long;
use Pod::Usage;
use File::Basename;

# Hardwired variables (can we keep them somewhere else?)

#my $cuffdiffPath = '/sw/rnaseq/cufflinks-2.0.2.Linux_x86_64';
my $cuffdiffPath = '/sw/opt/cufflinks-2.1.1';
my $samtoolsPath = '/sw/samtools-0.1.18';

my @conds = qw(WT Snf2);

my ($datapath, $outpath, $genome, $annotation, $conditions, $outfile);
my ($exrep, $inrep);  # unofficial options, for internal applications only
my ($help, $man);
my $threads = 1;

GetOptions(
   'datapath|d=s' => \$datapath,
   'outpath|p=s' => \$outpath,
   'outfile|o=s' => \$outfile,
   'genome|g=s' => \$genome,
   'annotation|a=s' => \$annotation,
   'conditions|c=s' => \$conditions,
   'exrep=s' => \$exrep, # unofficial
   'inrep=s' => \$inrep, # unofficial
   'threads=i' => \$threads,
   'man'       => \$man,
   'help|?'    => \$help,
);

pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;
die "Please supply a valid data path\n" unless $datapath && -d $datapath;
die "Please supply a valid output path\n" unless $outpath && -d $outpath;
die "Please supply a valid annotation file\n" unless $annotation && -e $annotation;
die "Please supply a valied genome file\n" unless $genome && -e $genome;
die "Need -outfile\n" unless $outfile;

@conds = split /,/, $conditions if defined $conditions;
my $labels = join ",", @conds;

# check condition directories
for my $cond (@conds)
  {die "$datapath/$cond does not exist\n" unless -d "$datapath/$cond"}

my %ex = map {$_, 1} split /,/, $exrep if $exrep;
my %in = map {$_, 1} split /,/, $inrep if $inrep;

my @samples = ();
for my $cond (@conds)
{
  my $dir = "$datapath/$cond";

  my @s = ();
  opendir DIR, $dir or die "Unable to open directory '$dir': $!\n";
  my @files = grep {/\.bam$/} readdir DIR;
  for my $file (@files)
  {
    my $base = basename($file, ".bam");
    my $rep = ($base =~ /rep(\d{2})/) ? $1 : '';
    $rep =~ s/^0+//;
    next if $exrep && $ex{$cond . $rep};
    next if $inrep && !$in{$cond . $rep};
    my $bam = "$dir/$file";
    push @s, $bam;
    #last if @s == 3;
  }
  closedir DIR;
  my $s = join ",", @s;      # all BAM files for one condition
  push @samples, $s
}
my $samples = join " ", @samples;    # two sets of BAM files

print <<EOF;

####################
  Running cuffdiff
(might take a while)
####################
EOF
#my $opt = "-o $outpath -u -b $genome -p 4 -L $labels";
my $opt = "-o $outpath -u -b $genome -p $threads -L $labels";
my $cmd = "$cuffdiffPath/cuffdiff $opt $annotation $samples";
#print "$cmd\n";die;
open R, ">$outpath/README" or die;
print R "Created with $cmd\n";
close R;
system($cmd) == 0 or die "Cuffdiff failed: $?";

print "\nCuffdiff successful (hopefully)\n";
my $exfile = "$outpath/gene_exp.diff";
open E, $exfile or die "Something went wrong: cufflinks DE file $exfile cannot be opened\n";
open F, ">$outfile" or die "Cannot open output file $outfile\n";
while(my $line = <E>)
{
  chomp $line;
  my @s = split /\t/, $line;
  my $gene = $s[1];
  my $fc = $s[9];
  my $p = $s[11];
  my $q = $s[12];
  my $sig = $s[13];
  my $e1 = $s[7];
  my $e2 = $s[8];
  print F "$gene\t$fc\t$p\t$q\t$sig\t$e1\t$e2\n";
}


=head1 SYNOPSIS

  run_cuffdiff.pl -d=/cluster/gjb_lab/cdr/GRNAseq/mapping/genome_biol_reps -p=cuffdiff_results -g=Scerevisiae68_ERCC92.fasta -a=Saccharomyces_cerevisiae.EF4.68.gtf -o=cuffdiff_expr.tsv

=head1 DESCRIPTION

This script is a wrapper around cuffdiff. It takes BAM files from provided directories, runs cuffdiff and converts its output to the required format.

=head1 OPTIONS

=over 5

=item B<-datapath|d>=C<path>

Path to the data directory of the experiment. Each sub-directory in this will be treated as a separate condition in the experiment (with the name of the condition matching the name of the directory), and each .bam file in each directory is a replicate for that condition.

=item B<-outpath|p>=C<path>

Path to the output directory, where cufflinks will write numerous output files.

=item B<-otufile|o>=C<filename>

The name (inc. path) of the output file from the wrapper.

=item B<-annotation|a>=C<filename>

Path to the .gff feature annotation file for the data. This file should match the feature you are counting the RNA-Seq expression for.

=item B<-genome|g>=C<filename>

Path to the FASTA file with reference genome. This is used by cuffdiff to run its bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates. You should have write access to the directory with this file as cufflinks will attempt to write an index file there. It might be a good idea to copy genome file (or a link) into your local directory.

=item B<-conditions|c>=C<string>

Comma-separated list of conditions. Default value is C<WT,Snf2>.

=item B<-threads>=I<number>

Number of threads to be used by cuffdiff.

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back



=cut