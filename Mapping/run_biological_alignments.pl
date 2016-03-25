#!/usr/bin/perl

=head1 NAME

run_biological_alignments.pl - given a set of combined biological replicate fastq files run tophat on them 

=cut

use strict;
use warnings;

use Getopt::Long qw(:config auto_version);
use Pod::Usage;
use File::Basename;
#use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

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
our $VERSION = '0.1';

GetOptions (
   'src-path=s' => \$inpath,
   'dest-path=s' => \$outpath,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid source path.') unless ($inpath && -d $inpath);
pod2usage(-msg => 'Please supply a valid destination path.') unless ($outpath && -d $outpath);

# check for required paths - both src and eest paths need to have 'WT' and 'Snf2' subdirs.
for my $d (qw/WT Snf2/) {
   die "ERROR - source dir doesn't have a '$d'/ subdir." unless ( -d "$inpath/$d");
   die "ERROR - destination dir doesn't have a '$d'/ subdir." unless ( -d "$outpath/$d");
}

# go through src dirs, find files and send tophat jobs to cluster.
print "Reading $inpath for fastq files...\n" if $VERBOSE;
my $n  =0;
for my $dir (qw/WT Snf2/) {
   opendir(my $dh, "$inpath/$dir") or die "ERROR - unable to open directory '$inpath/$dir': $!\nDied";
   while (my $entry = readdir($dh)) {
      next if ($entry =~ /^\./);  # skip hidden/special files
      my $base = basename($entry,".fastq.gz");
      my $rep = '';
      if ($base =~ /rep(\d{2})/) {
         $rep = $1;
      }
      print "Submitting Tophat run on $entry\n" if $VERBOSE;
      my $CMD = "source /gridware/sge/default/common/settings.sh && qsub -q 64bit-pri.q -pe smp 4 -v PATH=$bowtiePath:$tophatPath:$samtoolsPath:\$PATH -N TH${dir}_r$rep -M c.cole\@dundee.ac.uk -m a -cwd -b y -R y tophat --num-threads 4 --output-dir $outpath/$dir/${base}_tophat2.0.5 $tophatParams $genome $inpath/$dir/$entry";
      print "DEBUG CMD: $CMD\n" if $DEBUG;
      system("$CMD") == 0 or die "ERROR - system() failed for '$CMD'\nDied";
      ++$n;
   }
   closedir($dh);
}
print "Submitted $n jobs to the cluster\n";

=head1 SYNOPSIS

run_biological_alignments.pl --src-path <path> --dest-path <path> [--verbose|--no-verbose] [--version] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

Use this to run tophat over the set of fastq files where technical replicates have been combined into biological replicates.

For this experiment it means 48 replicates each of WT and Snfs2 samples. 

=head1 OPTIONS

=over 5

=item B<--src-path>

Path containing 'WT' and 'Snf2' subdirs with fastq files.

=item B<--dest-path>

Path containing 'WT' and 'Snf2' subdirs for saving tophat output.

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

Chris Cole <christian@cole.name>

=head1 COPYRIGHT

Copyright 2012, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut