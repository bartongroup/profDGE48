#!/usr/bin/perl

=head1 NAME

run_alignments.pl - take raw Fastq files and align against genome sequence via submission to cluster

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $inpath;
my $outpath;
my $genome = '/db/bowtie2/Scerevisiae68_ERCC92';
my $tophatParams = '--max-intron-length 1000 --min-intron-length 10 --microexon-search --b2-very-sensitive --max-multihits 1';
my $tophatPath = '/sw/rnaseq/tophat-2.0.5.Linux_x86_64';
my $bowtiePath = '/sw/rnaseq/bowtie2-2.0.0-beta7';
my $samtoolsPath = '/sw/samtools-0.1.18';
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
my $version = '0.2';

$| =1;
GetOptions (
   'src=s'     => \$inpath,
   'dest=s'    => \$outpath,
   'genome=s'  => \$genome,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid input path.') unless ($inpath && -d $inpath);
pod2usage(-msg => 'Please supply a valid output path.') unless ($outpath && -d $outpath);

die "ERROR - output path doesn't have a 'WT' subdir\n" unless (-d "$outpath/WT");
die "ERROR - output path doesn't have a 'Snf2' subdir\n" unless (-d "$outpath/Snf2");

print "\nSrc dir: $inpath\n";
print "Dest dir: $outpath\n";
print "Genome: $genome\n";
print "Samtools path: $samtoolsPath\n";
printf "Samtools version: %s\n", samtoolsVer($samtoolsPath);
print "Tophat path: $tophatPath\n";
printf "Tophat version: %s\n", tophatVer($tophatPath);
print "Bowtie path: $bowtiePath\n";
printf "Bowtie version: %s\n", bowtieVer($bowtiePath);
print "Tophat settings: $tophatParams\n\n";


## assumes there are dirs named '1' to '96' in source path and in each dir there are several fastq files,
## and that 1-48 are Snf2 mutants and 49-96 are WTs.

# do Snf2 files first
print "Searching for and submitting Snf2 Fastq files...\n" if $VERBOSE;
my $tot = 0;
my $cond = 'Snf2';
foreach my $n (1..48) {
   opendir(my $dh, "$inpath/$n") or die "ERROR - unable to open directory '$inpath/$n': $!\n";
   my $i = 1;
   while(my $entry = readdir($dh)) {
      next if ($entry =~ /^\./); # skip hidden/special files
      if ($entry =~ /fastq\.gz$/) {
         my $CMD = sprintf("source /gridware/sge/default/common/settings.sh && qsub -q 64bit-pri.q -pe smp 4 -v PATH=$bowtiePath:$tophatPath:$samtoolsPath:\$PATH -N TH${cond}_r%02d_l%d -M c.cole\@dundee.ac.uk -m a -cwd -b y -R y tophat --num-threads 4 --output-dir $outpath/$cond/rep%02d_lane%d $tophatParams $genome $inpath/$n/$entry", $n, $i, $n, $i);
         print "DEBUG CMD: $CMD\n" if $DEBUG;
         system($CMD) == 0 or die "ERROR - system() failed for '$CMD'\nDied";
         #print "$CMD\n";
         ++$tot;
      }
      ++$i;
   }
   closedir($dh);
}

# then WT files...
print "\nSearching for and submitting WT Fastq files...\n" if $VERBOSE;
$cond = 'WT';
foreach my $n (49..96) {
   opendir(my $dh, "$inpath/$n") or die "ERROR - unable to open directory '$inpath/$n': $!\n";
   my $i = 1;
   while(my $entry = readdir($dh)) {
      next if ($entry =~ /^\./); # skip hidden/special files
      if ($entry =~ /fastq\.gz$/) {
         my $CMD = sprintf("source /gridware/sge/default/common/settings.sh && qsub -q 64bit-pri.q -pe smp 4 -v PATH=$bowtiePath:$tophatPath:$samtoolsPath:\$PATH -N TH${cond}_r%02d_l%d -M c.cole\@dundee.ac.uk -m a -cwd -b y -R y tophat --num-threads 4 --output-dir $outpath/$cond/rep%02d_lane%d $tophatParams $genome $inpath/$n/$entry", $n-48, $i, $n-48, $i);
         print "DEBUG CMD: $CMD\n" if $DEBUG;
         system($CMD) == 0 or die "ERROR - system() failed for '$CMD'\nDied";
         ++$tot;
      }
      ++$i;
   }
   closedir($dh);
}
print "Submitted $tot files for Tophat alignment...\nCheck with qstat to see when jobs finish\n" if $VERBOSE;

# quick sub to extract tophat version string
sub tophatVer {
   my $path = shift;
   
   my $out = `$path/tophat --version` or die "ERROR - unable to execute '$path/tophat --version': $!\n";
   if ($out =~ /v([\d\.]+)/) {
      return($1);
   } else {
      die "ERROR - unable to capture version information for Tophat\n";
   }
}

#quick sub to extract bowtie version string
sub bowtieVer {
   my $path = shift;
   
   my $out = `$path/bowtie2 --version` or die "ERROR - unable to execute '$path/tophat --version': $!\n";
   if ($out =~ /version ([\d\.abet\-]+)/) {
      return($1);
   } else {
      die "ERROR - unable to capture version information for Tophat\n";
   }
   
}

#quick sub to extract samtools version string
sub samtoolsVer {
   my $path = shift;
   
   my $out = `$path/samtools 2>&1` or die "ERROR - unable to execute '$path/tophat --version': $!\n";
   if ($out =~ /Version: ([\d\.]+)/) {
      return($1);
   } else {
      die "ERROR - unable to capture version information for Tophat\n";
   }
   
}

=head1 SYNOPSIS

run_alignments.pl --src <dir> --dest <dir> [--genome <pathprefix>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

This scripts takes all 96 replicates from the yeast snf2 mutant experiment and aligns the raw Fastq reads to the genome with Tophat.

The i<--src> directory must have 96 subdirs numbered 1-96, where 1-48 are the snf2 mutants and 49-96 are the wild-types.

The i<--dest> directory must have the Snf2/ and WT/ subdirs where the Tophat runs will be saved.

=head1 OPTIONS

=over 5

=item B<--src>

Source path of subdirs with fastq files 

=item B<--dest>

Destination path for output.

=item B<--genome>

Path and to bowtie2 genome index [default: /db/bowtie2/Scerevisiae68_ERCC92]

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

Chris Cole <christian@cole.name>

=head1 COPYRIGHT

Copyright 2012, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut