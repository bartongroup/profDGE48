package GetOpt;

=head1 NAME

GetOpt 

=head1 SYNOPSIS

  use GetOpt;
  
=head1 DESCRIPTION
  

=cut

require Exporter;

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use GRNASeq;

our @ISA = ("Exporter");
our @EXPORT = qw(
  $psfile
  $cond
  $rep
  $nonzero
  $norm
  $type
  $exclude
  $include
  $clip
  $genes
  $genlist
  $gff
  $clean
  $spclean
  $randrep
  $randrep2
  $multicor
  $help
  $man
  GetMyOptions
);


our $psfile;
our $cond = 'WT';
our $rep;
our $nonzero = 1;
our $norm = 'deseq';
our $type = 'raw';
our $exclude;
our $include;
our $clip;
our $genes;
our $genlist;
our $clean;
our $spclean;
our $randrep;
our $randrep2;
our $multicor = 'bh';  # Benjamini-Hochberg
our $gff;
our ($help, $man);


sub GetMyOptions
{
  my %opt = @_;

  my $stat = GetOptions(
    'psfile=s' => \$psfile,
    'cond=s' => \$cond,
    'rep=i' => \$rep,
    'nonzero=i' => \$nonzero,
    'norm=s' => \$norm,
    'type=s' => \$type,
    'exclude=s' => \$exclude,
    'include=s' => \$include,
    'clip=i' => \$clip,
    'genes=s' => \$genes,
    'genlist=s' => \$genlist,
    'randrep=i' => \$randrep,
    'randrep2=i' => \$randrep2,
    'multicor=s' => \$multicor,
    'countsdir=s' => \$countsdir,
    'outlierdir=s' => \$outlierdir,
    'pileupdir=s' => \$pileupdir,
    'gff=s' => \$gff,
    clean => \$clean,
    spclean => \$spclean,
    help => \$help,
    man => \$man,
    %opt
  );
  die "Failed getting options. Check your command-line options and try again.\n" unless $stat;
  pod2usage(-verbose => 2) if $man;
  pod2usage(-verbose => 1) if $help;
}

=head1 GLOBAL OPTIONS

These are options in the GetOpt module, shared by many scripts.

=over 4

=item B<-psfile>=I<pathname>

Output postscript file.

=item B<-cond>=I<string>

Condition to use. In some scripts analysis can be restricted to one condition.

=item B<-rep>=I<number>

Replicate number.

=item B<-nonzero>=I<number>

Selects only data rows (genes) with at least I<n> non-zero replicates.

=item B<-norm>=I<string>

Expression data normalization. Input expression files contain raw counts. These are normalized on the fly, using one of the following methods:

  none - no normalization
  deseq - DESeq normalization (default)
  tmm - trimmed mean of M values normalization
  fpkm - approximate FPKM normalization (not identical to cuffdiff!)
  totcount - total count
  totcountm - total count with normalization factors scaled to the mean of 1

=item B<-type>=I<string>

Type of raw data, as encoded in the counts file name (e.g. WT_raw.tsv is 'raw'). Could be 'lev' for levelled data, or some other, alternative prescription.

=item B<-exclude>=I<string>

A list of replicates to exclude from analysis. The format of this string can be, e.g.
   
   "1,3,5:" to exclude replicates 1, 3 and 5 from the first condition
   ":1,2,3 "to exclude replicates 1, 2 and 3 from the second condition
   "1:3,4" to exclude replicates 1 from the first condition and 3 and 4 from the second contition  

=item B<-include>=I<string>

A list of replicates to include to analysis. Format as above.

=item B<-clip>=I<number>

Clip C<n> top and n bottom replicates.

=item B<-genes>=I<string>

Comma-delimeted list of genes.

=item B<-genlist>=I<pathname>

A file with a list of genes (first column).

=item B<-randrep, -randrep2>

Internal options.

=item B<-clean>

A switch indicating that only clean replicates should be used. Clean replicates are defined in C<GRNASeq.pm> module.

=item B<-multicor>=[none|BH|HB]

Multiple test correction to apply to tests/tools that return uncorrected raw p-values. The default value is C<BH>.
 
=item B<-countsdir>=I<path>
 
Path to the directory containting gene read count data. There should be two files, one for each conditions, with names C<cond>_raw.tsv, created by C<combine_replicates.pl>. If not specified, the value from <defs.dat> will be used.

 =item B<-outlierdir>=I<path>
 
 A directory to store outlier data. If not specified, the value from <defs.dat> will be used.
 
 =item B<-pileupdir>=I<path>
 
 A directory to store pileup data. If not specified, the value from <defs.dat> will be used.
 
 =item B<-gff>=I<pathname>
 
 A GTF file to be used by many scripts. If not specified, the value from <defs.dat> will be used.
 
 =item B<-clean>

A switch indicating that only clean replicates should be used. Clean replicates are defined in C<GRNASeq.pm> module.
 
 =item B<-spclean>

A switch indicating that only spikein-clean replicates should be used. Spkein-clean replicates are defined in C<GRNASeq.pm> module.
 
 

=back

=cut

1;
