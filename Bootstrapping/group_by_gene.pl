#!/usr/bin/perl

=head1 NAME

group_by_gene.pl - given a BAM file and a GTF file group reads by gene

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Temp;

my $bam;
my $gtf;
my $out = 'genes.out';
my $method = 'union';
my $feature = 'gene_id';
my $useLocal = 1;
my $VERBOSE = 1;
my $DEBUG = 0;
my $JAVA = '/sw/java/jdk1.6.0_32/bin/java';
my $picardPath = '/homes/ccole/bin';
my $SAMTOOLS = '/sw/samtools-0.1.18/samtools';
my $HTSEQ = '/sw/bin/htseq-count';
my $tool = 'samtools';
my $help;
my $man;
my $version = "1.8.2";

# force unbuffering of print streams
$|=1;

GetOptions (
   'bam=s'     => \$bam,
   'gtf=s'     => \$gtf,
   'out=s'     => \$out,
   'method=s'  => \$method,
   'feature=s' => \$feature,
   'samtools=s' => \$SAMTOOLS,
   'tool=s'    => \$tool,
   'htseq=s'   => \$HTSEQ,
   'local!'    => \$useLocal,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
   'version'   => sub { printVers(); exit }
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid BAM filename.') unless ($bam && -s $bam);
pod2usage(-msg => 'Please supply a valid GTF filename.') unless ($gtf && -s $gtf);
pod2usage(-msg => 'Please supply a valid method [union|istrict|inoempty].') unless ($method eq 'union' or $method eq 'istrict' or $method eq 'inoempty');
pod2usage(-msg => 'Please supply a valid feature [gene_id|transcript_id].') unless ($feature eq 'gene_id' or $feature eq 'transcript_id');

# print out version and path info for this script, samtools and htseq-count
printVers();

printf "HOST: %s\n", $ENV{HOSTNAME} if $DEBUG;

# strip path from BAM filename
my $filename = basename($bam);

## First convert BAM into SAM file, but sorted by readID (as required by htseq-count)
## using 'samtools'.
## Must use the local filesystem as the SAM file can be large
print "Creating sorted SAM file for 'htseq-count...\n" if $VERBOSE;

# Determine 'local' filesystem path and create a temporary file prefix there for use during sorting
my $local = './';
if ($useLocal) {
   if (defined($ENV{'TMPDIR'})) {
      $local = $ENV{'TMPDIR'};
   } else {
      warn "Warning - using './' as temporary path as TMPDIR not set\n";
   }
}
my $tmp = File::Temp->new( DIR => $local );
my $tmpFile = $tmp->filename();
print "Temporary sorting file prefix: $tmpFile\n" if $DEBUG;

# create temporary SAM file on 'local' filesystem
my $sam = "$local/$filename.sam";
print "Temporary output file: $sam\n" if $DEBUG;

# run samtools... 
print "Sorting BAM file...\n" if $VERBOSE;
my $CMD = '';
if ($tool eq 'samtools') {
   $CMD = "$SAMTOOLS sort -n $bam $tmpFile";
   print "executing $CMD \n" if $DEBUG;
   system($CMD) == 0 or die "ERROR - system() failed for command: $CMD\nDied";
   die "ERROR - no output or output is empty from 'samtools sort $tmpFile.bam'\n" unless ("$tmpFile.bam" && -s "$tmpFile.bam");
   
   print "Converting BAM -> SAM...\n" if $VERBOSE;
   $CMD = "$SAMTOOLS view -h $tmpFile.bam > $sam";
   print "executing $CMD \n" if $DEBUG;
   system($CMD) == 0 or die "ERROR - system() failed for command: $CMD\nDied";
   die "ERROR - no output or output is empty from 'samtools view'\n" unless ($sam && -s $sam);
} elsif ($tool eq 'picard') {
   $CMD = "$JAVA -Xms1000m -Xmx1500m -jar $picardPath/SortSam.jar INPUT=$bam OUTPUT=$sam SORT_ORDER=queryname VALIDATION_STRINGENCY=LENIENT";
   print "executing $CMD \n" if $DEBUG;
   system($CMD) == 0 or die "ERROR - system() failed for command: $CMD\nDied";
   die "ERROR - no output or output is empty from 'SortSam.jar $bam'\n" unless ("$sam" && -s "$sam");
} else {
   die "ERROR - unable to sort with '$tool' tool. Unknown tool. Please specify 'samtools' or 'picard'\n"
}

## Then use SAM file to group reads to genes using 'htseq-count'
print "Running htseq-count to generate gene read counts...\n" if $VERBOSE;
if ($method eq 'istrict') {
   $CMD = "$HTSEQ --mode 'intersection-strict' --stranded 'no' --type 'exon' --idattr '$feature' $sam $gtf -q > $out";
} elsif ($method eq 'inoempty') {
   $CMD = "$HTSEQ --mode 'intersection-nonempty' --stranded 'no' --type 'exon' --idattr '$feature' $sam $gtf -q > $out";
} else {
   $CMD = "$HTSEQ --mode 'union' --stranded 'no' --type 'exon' --idattr '$feature' $sam $gtf -q > $out";
}
print "executing $CMD \n";
system($CMD) == 0 or die "ERROR - system() failed for command: $CMD\nDied";
die "ERROR - no output or output is empty from '$HTSEQ'\n" unless ($out && -s $out);
unlink($sam, "$tmpFile.bam");
print "Done!\n";

# get the version string for a given program name
sub getVersion {
   my $prog = shift;
   
   my $version = '';
   if ($prog =~ /samtools$/) {
      my $out = `$prog 2>&1`;
      if ($out =~ /Version: (.*)/) {
         $version = $1;
      }
   } elsif ($prog =~ /htseq-count$/) {
      my $out = `$prog 2>&1`;
      if ($out =~ /version (.*)\./) {
         $version = $1;
      }
   } elsif ($prog =~ /java$/) {
      my $out = `$prog -Xms150m -version 2>&1`;
      if ($out =~ /version \"(.*)\"/) {
         $version = $1;
      }
   } elsif ($prog =~ /jar$/) {
      my $out = `$JAVA -Xms150m -jar $prog 2>&1 --version`;
      if ($out =~ /(.*)/) {
         $version = $1;
      }
   } else {
      die "ERROR - program '$prog' not known. Can't extract version number from it\nFix the code, man!\nDied";
   }
   
   if ($version) {
      return($version)
   } else {
      warn "Warning - unable to get version info for '$prog'\n";
      return(0)
   }
}

sub printVers {
   
   print "$0 version: $version\n";
   die "ERROR - htseq-count path '$HTSEQ', non-existent or not executable\n" unless (-e $HTSEQ && -x $HTSEQ);
	printf "$HTSEQ version: %s\n", getVersion($HTSEQ);
   if ($tool eq 'samtools') {
      die "ERROR - samtools path '$SAMTOOLS', non-existent or not executable\n" unless (-e $SAMTOOLS && -x $SAMTOOLS);
      printf "$SAMTOOLS version: %s\n", getVersion($SAMTOOLS);
   } else {
      die "ERROR - java path '$JAVA', non-existent or not executable\n" unless (-e $JAVA && -x $JAVA);
      die "ERROR - SortSam path '$picardPath/SortSam.jar', is non-existent\n" unless (-e "$picardPath/SortSam.jar");
      printf "$JAVA version: %s\n", getVersion($JAVA);
      printf "$picardPath/SortSam.jar version: %s\n", getVersion("$picardPath/SortSam.jar");
   }
   return(1);

}

=head1 SYNOPSIS

group_by_gene.pl --bam <file> --gtf <file> [--tool <samtools|picard>] [--method <union|istrict|inoempty>] [--feature <gene_id|transcript_id>] [--out <file>] [--samtools <path>] [--htseq <path>] [--local|--no-local] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help] [--version]

=head1 DESCRIPTION

Essentially, this script is a wrapper for F<htseq-count>, but makes it useable with BAM files.

F<htseq-count> only works with SAM files sorted by readID whereas most other utilities require BAM files sorted by chromosomal position. This script takes a BAM file converts it into an appropriately sorted SAM file and runs F<htseq-count> on it.

As SAM files can be pretty large and are only required temporarily here, the SAM file is stored on the TMPDIR path and deleted. Thereby, saving precious GPFS disk space. This behaviour can be controlled with the B<--local|--no-local> switch. 

The script is optimised for the cluster so use the B<ram> and B<local_free> I<SGE> options to ensure the jobs don't run out of local disk space nor RAM:

   qsub -q 64bit-pri.q -l ram=4G -l local_free=70G group_by_gene.pl --bam my.bam --gtf my.gtf --out my.out

=head1 REQUIREMENTS

'htseq-count' and 'samtools' are required and must be in the default paths or provided with the --samtools and --htseq arguments.

=head1 OPTIONS

=over 5

=item B<--bam>

Input BAM file.

=item B<--gtf>

Input GTF file.

=item B<--tool>

Select which tool to use for sorting (samtools or picard). [default: samtools]

=item B<--method>

Specify which method to use for aggregating reads to genes. [default: union]

=item B<--feature>

Specify which feature to aggregate read on. [default: gene_id]

=item B<--out>

Output filename. [default: gene.out]

=item B<--samtools>

Give path to samtools executable [default: /sw/samtools-0.1.18/samtools]

=item B<--htseq>

Give path to samtools executable [default: /sw/bin/htseq-count]

=item B<--local|--no-local>

Toggle whether to use /local for tmp files, otherwise use current dir. [default --local]

=item B<--verbose|--no-verbose>

Toggle verbosity. [default:none]

=item B<--debug|--no-debug>

Toggle debugging output. [default:none]

=item B<--version>

Print version info and exit

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back

=head1 AUTHOR

Chris Cole <christian@cole.name>

=cut