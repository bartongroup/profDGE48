#!/sw/bin/perl

=head1 NAME

B<pileup_var.pl>

=head1 DESCRIPTION

Create files with pileup statistic for each gene: one file with mean and one with standard deviation. There is one gene per line, and each line contains gene name and 48 numbers for all replicates.

Requires pileup files for each chromosome, created with C<build_pileup.pl>. 

=cut


use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use PDL;
use PDL::NiceSlice;
use PDL::Graphics::PLplot;

use GRNASeq;

use PLgraphs;
use Tools;

Stamp();
$| = 1;

my $gff = 'Saccharomyces_cerevisiae.EF4.68.gtf';
my $outdir;
my ($help, $man);
GetOptions(
  'pileupdir=s' => \$pileupdir,
  'outdir=s' => \$outdir,
  'gff=s' => \$gff,
  help => \$help,
  man => \$man
);
pod2usage(-verbose => 2) if $man;
pod2usage(-verbose => 1) if $help;

die "Provide valid GFF file\n" unless defined $gff && -e $gff;
die "Provide valid pileup directory\n" unless defined $pileupdir && -d $pileupdir;

$outdir ||= $pileupdir;

#my $gfile = CountFile('yeast', 'WT', undef, 'raw');
#my @genes = ReadGeneList($gfile);
my ($gendat, $genidx) = ReadGFF($gff);
my @genes = keys %$gendat;

# read gene data, arrange by chromosomes
my %genes = ();
for my $gene (@genes)
{
  next unless defined $gendat->{$gene};
  my ($chr, $start, $end, $strand, $gene_name) = GeneStructure($gene, $gendat, 1);
  $genes{$chr}{$gene} = [$start, $end, $strand, $gene_name];
  
  #print "$gene $chr\n";
}

for my $cond (@conds)
{
  my $e = ($gff =~ /ERCC/) ? '_ERCC' : '';
  my $mfile = "$outdir/${cond}${e}_pileupvar_mean.tsv";
  my $sfile = "$outdir/${cond}${e}_pileupvar_sd.tsv";
  #my $chifile = "$outdir/${cond}${e}_pileupvar_chisq.tsv";
  open M, ">$mfile" or die;
  open S, ">$sfile" or die;
  #open X, ">$chifile" or die;
  for my $chr (keys %genes)
  {
    my $pfile = "${pileupdir}/${cond}_chr${chr}_pileup.tsv";
    die "Cannot find file $pfile\n" unless -e $pfile;

    local *F;
    open F, $pfile or die;
    my $line = <F>;
    my @s = split " ", $line;
    my $ncol = @s;
    
    print "Reading pileup file for $cond/$chr...";
    my ($pos, @d) = rcols $pfile, (0 .. $ncol-1); 
    my $d = pdl(@d)->transpose();
    my ($nrep, $npos) = dims($d); 
    print " done\n";
    
    #my ($mm) = statsover($d); my ($mmm) = stats($mm); print "### $mmm\n"; 
    
    my @chrgenes = (keys %{$genes{$chr}});
    my $N = @chrgenes;
    my $cnt = 0;
    print "Processing genes...       ";
    for my $gene (keys %{$genes{$chr}})
    {
      Percent(($cnt++)/$N);
      print M "$gene";
      print S "$gene";
      my ($start, $end) = @{$genes{$chr}{$gene}};
      for my $rep (1 .. $nrep)
      {
        my $gdat = $d($rep-1,$start-1:$end-1;-);
        #print "$cond $chr $gene $rep $gdat\n";
        my ($m, $s) = stats($gdat);
        my $n = nelem($gdat);
        #my $chisq = ($m > 0) ? $s**2 / $m : 0;   # reduced chi-squared
        #my $cv = ($m > 0) ? $s/$m : 0;
        printf M "\t%.4g", $m;
        printf S "\t%.4g", $s;
        #printf X "\t%.4g", $chisq;
      }
      print M "\n";
      print S "\n";
      #print X "\n";
    }
    print "\b\b\b\b\b\b\b  done    \n";
  }
}


sub ReadPileup
{
  my ($cond, $chr, $pos1, $pos2, %opt) = @_;
  
  my $pfile = "${pileupdir}/${cond}_chr${chr}_pileup.tsv";
  die "Cannot find file $pfile\n" unless -e $pfile;
  
  my $nf = rcols "${cond}_$opt{norm}_normfac.txt" if defined $opt{norm};
  
  my ($cl1, $cl2) = split /:/, $CleanExclude;
  my %cl = (WT => $cl1, Snf2 => $cl2);  
  $opt{exclude} = $cl{$cond} if $opt{clean};
  my $ex = pdl(split /,/, $opt{exclude}) - 1 if defined $opt{exclude};
  
  local *F;
  open F, $pfile or die;
  my $line = <F>;
  my @s = split " ", $line;
  my $ncol = @s;
  my $lines = $pos1-1 . ":" . $pos2;
  
  print "Reading pileup file...";
  my ($pos, @d) = rcols $pfile, (0 .. $ncol-1), {LINES=>$lines}; 
  print " done\n";
    
  my $d = pdl(@d)->transpose();
  my ($nrep, $npos) = dims($d); 
    
  $d /= $nf if defined $nf;
  
  if(defined $ex)
  {  
    my $all = sequence($nrep);
    my $xsel = setops($all, 'XOR', $ex);
    $d = $d($xsel,);
  }
  
  my $show;
  if(defined $opt{show})
  {
    $show = $d($opt{show}-1,;-)
  }
  
  my ($m, $s) = statsover($d);
  
  return $pos, $m, $s, $show
}



=head1 SYNOPSIS

  pileup_var.pl

=head1 OPTIONS

=over 5

=item B<-pileupdir>=I<pathname>

The directory with pileups (created with C<build_pileup.pl>). If not defined, the default value from C<defs.dat> file will be used.


=item B<-outdir>=I<pathname>

The directory to save results. If not specified, it will the pileup directory.

=item B<-gff>=I<pathname>

The GTF file to be used. The default value is 'Saccharomyces_cerevisiae.EF4.68.gtf'. 

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back

=cut
