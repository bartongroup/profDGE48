#!/usr/bin/perl

=head1 NAME

B<buil_pileup.pl>

=head1 DESCRIPTION

Build aggregated pileup from all BAM files. 

=cut

use strict;
use warnings;

$|=1;

use Getopt::Long qw(:config auto_version);
use Pod::Usage;
use File::Basename;
use File::Temp qw(tempfile);
use PDL;
use PDL::NiceSlice;

my $samtoolsPath = '/sw/samtools-0.1.18';

#my %c = ReadChromosomeLength('q1.bam');
#for my $c (keys %c)
#  {print $c, "->", $c{$c}, "\n"}
#die;

#my %p = GetPileup('q1.bam');
#my ($p, $c) = @{$p{I}};
#wcols $p, $c;die;

my @conds = qw(WT Snf2);
#my @chr = qw(I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI);

my ($datapath, $outpath, $genome, $conditions);
my ($help, $man);
my $spikeins;

GetOptions(
   'datapath|d=s' => \$datapath,
   'outpath|o=s' => \$outpath,
   'conditions|c=s' => \$conditions,
   'genome|g=s' => \$genome,
   'man'       => \$man,
   'help|?'    => \$help,
   spikeins => \$spikeins
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;
die "Please supply a valid data path\n" unless $datapath && -d $datapath;
die "Please supply a valid output path\n" unless $outpath && -d $outpath;
#die "Please supply a valied genome file\n" unless $genome && -e $genome;

@conds = split /,/, $conditions if defined $conditions;
my $labels = join ",", @conds;

# check condition directories
for my $cond (@conds)
  {die "$datapath/$cond does not exist\n" unless -d "$datapath/$cond"}

my %chrlen;
my $first = 1;
for my $cond (@conds)
{
  my $dir = "$datapath/$cond";
  my %pileup = ();
  opendir DIR, $dir or die "Unable to open directory '$dir': $!\n";
  my @files = grep {/\.bam$/} readdir DIR;
  my $n = 0;
  for my $file (@files)
  {
    my $bamfile = "$dir/$file";
    if($first)
    {
      %chrlen = ReadChromosomeLengths($bamfile);
      $first = 0;
    }
    #print join(", ", keys %chrlen), "\n";die;
    
    print "Reading $file...\n";
    my %p = GetPileup($bamfile, \%chrlen);
    print "\n";
    for my $chr (keys %p)
      {push @{$pileup{$chr}}, $p{$chr}}
    $n++;
    #last if $n==2;
  }
  closedir DIR;
  
  print "Processing...";
  for my $chr (keys %chrlen)
  {
    # 2D piddle: cols-replicates, rows-positions
    my $d = transpose(pdl(@{$pileup{$chr}}));
    my ($nrep, $npos) = dims($d);
    my $pos = sequence($npos) + 1;
    $d = append(transpose($pos), $d);   # add column of positions
    
    my $outfile = "$outpath/${cond}_chr${chr}_pileup.tsv";
    my @s = ();
    push @s, $d($_,;-) for (0 .. $nrep);  # list of columns
    wcols @s, $outfile;                   # cannot wcols 2D piddle 
  }
  print "done\n";
}


######################################################

# this works only for one particular genome!
sub ReadChromosomeLengths
{
  my ($bam) = @_;
  
  die "Cannot open BAM file $bam\n" unless -e $bam;

  my ($F, $file) = tempfile('bamXXXXX', DIR=>'.', SUFFIX=>'.txt', UNLINK=>1);
  close $F;
  `$samtoolsPath/samtools view -H $bam > $file`;
  die "$?\n" if $?;
  
  my %h = ();
  local *F;
  open F, $file or die;
  while(my $line = <F>)
  {
    chomp $line;
    next if (!$spikeins && $line =~ /ERCC/) || ($spikeins && $line !~ /ERCC/);   # spikeins
    if($line =~ /^\@SQ/)
    {
      my @s = split " ", $line;
      my ($sn, $ln);
      for my $s (@s)
      {
        if($s =~ /SN:(.*)/)
          {$sn = $1}
        if($s =~ /LN:(.*)/)
          {$ln = $1}
      }
      die "Unrecognized SQ line in $bam:\n$line\n" unless defined $sn && defined $ln;
      $h{$sn} = $ln
    }
  }
  return %h
}

sub GetPileup
{
  my ($bam, $chrlen) = @_;
  
  # create safe temp file ('tempfile' opens the file, so need to close it)
  my $dir = '.';
  my ($F, $file) = tempfile('pileupXXXXX', DIR=>$dir, SUFFIX=>'.txt', UNLINK=>1);
  close $F;
  
  # prepare zero tables for each chromosome
  my %p = ();
  for my $chr (keys %$chrlen)
    {$p{$chr} = zeroes($chrlen->{$chr})}
    
  #system("$samtoolsPath/samtools mpileup -D $bam > $file") == 0 or die;
  `$samtoolsPath/samtools mpileup -D $bam > $file`;
  die if $?;
  local *F;
  open F, $file or die;
  while(my $line = <F>)
  {
    chomp $line;
    my ($chr, $pos, $base, $cnt) = split " ", $line;
    next unless defined $chrlen->{$chr};   # skip spikeins
    if($pos <= $chrlen->{$chr})
      {$p{$chr}->($pos - 1) .= $cnt}
    else
      {die "Position $pos outside chromosome $chr\n"}
  }
  close F;
  unlink $file;
  
  return %p
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

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back



=cut