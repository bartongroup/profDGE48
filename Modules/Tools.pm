package Tools;
require Exporter;
require PDL;
require PDL::NiceSlice;
require Term::Complete;
require File::Copy;
require File::Temp;
require Devel::Size;


use Term::Complete;
use File::Copy;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use Devel::Size qw(total_size);
use File::Temp qw(tempfile mktemp);

our $VERSION = 1.00;

our @ISA = ("Exporter");
our @EXPORT = qw(
 call warn ToBe
 Percent
 Percent1
 Percent2
 ReadASCII
 ReadFileList
 ReadLine
 ReadName
 ReadYesNo
 Retrieve
 Rmmin
 Store
 Stamp
 sStamp
 WriteASCII
 ReadFastaSequences
 WriteFastaSequences
 SmoothMovingMean
 wc
 commify
);

################################################################  

sub call
#
# Call to a system task
#
# @log = call($arguments, $status [, $ignore])
#
# Call a system command (e.g. an ftool) and process errors. Returns $status = $?,
# i.e. 0 if there is no error.
#
# If $ignore == 1 then warnings are ignored (i.e. $? = 0, but there is something
# sent to STDERR)
#
{
  my $command = $_[0];
  my $ignore = (defined $_[2]) ? $_[2] : 0;
  
  my $tmp = defined $ENV{TMP} ? $ENV{TMP} : '/tmp'; 
  my $errfile = "$tmp/" . mktemp("errXXXXXX"); 
  unlink $errfile;
  
  my @result = `$command 2> $errfile`;
  
  $_[1] = $?;
  if($?)
  {
    print "\n### This command went wrong:\n";
    print "$command\n\n";
  }
  my $errmsg = '';
  if(open ERR, $errfile)
  {
    while(my $line = <ERR>)
    {
      chomp $line;
      $errmsg .= "$line\n" if $line !~ /^\s*$/;
    }
  }
  if($errmsg && !($ignore && !$?))
  {
    print "\n### Error output:\n";
    print "$errmsg\n";
  }
  unlink $errfile;
  return @result;
}

sub call_old {

# @log = &call($args,"warn","Problem running:\n $args","logfile.log")

## based on runcom, but keeps the output from the ftool even in the case
## of an error;  give command, error routine, error string, error log;
## returns array:  1st value is 0 if good or $? or 1 if system or ftool error; 
## then command; then error messages sent by ftool or system to 
## STDERR; then any other output from tool;  retains suppress option.

    my($command) = $_[0];
    my($ERRFIL,@result,$CURHANDLE,$errormsg);
    my $suppress;
    
    if($command =~ s/^-// ) {$suppress = 1;}
    
    $ERRFIL = "/tmp/error$$";
    unlink $ERRFIL;
    $command .= " 2> $ERRFIL";
    @result = `$command`;

    if ( !open(ERRFIL, $ERRFIL) ) {  ## if no ERRFIL created, return with $? status
        unlink $ERRFIL;
    }
    else {
LOP:  while($_ = <ERRFIL>)
             {
	       last LOP if !(/(ld\.so|^\s*\n$)/);
             }
      if ( ! $_ ) {  ## if ERRFIL created but empty, return
	  unlink $ERRFIL;
      }
      else {  ##  there is an error message
	  unless ($suppress) {
	    $errormsg = $_;
            unlink $ERRFIL;
            }
      }  ## else there's an error message
    } ##  else there's an ERRFIL created

    if ($? || ($errormsg))        ## call err_routine if either $? or $errormsg
    { 
      $CURHANDLE = select(STDERR);
      print "\n-------------------------------------------------------------\n";
      print "Something went wrong:\n";
      print "$command\n";
      print "-------------------------------------------------------------\n";
      select($CURHANDLE);

      if($errormsg)
      {
        push @result, "Problem running:\n $command:\n\n$errormsg\n";
        $_[1] = 1;
      }
      else
      {
        push @result, "Problem running:\n $command\n";
        $_[1] = $?;
      }
    } ## if any errors

    else 
    {
      push @result, "$command\n";
      $_[1] = $?;
    }
    return @result;

} ## sub call

sub warn {

    my $string = $_[0];
    my $log = $_[1];
    my $CURHANDLE = select(STDERR);
    print "-------------------------------------------------------------\n";
    print "$string\n";
    print "See full log in $log\n" if $log;
    print "-------------------------------------------------------------\n";
    select($CURHANDLE);
}

sub ToBe
{
  my $stat = $_[0];
  my $log = $_[1];

  if($stat)
  {
    die("I'm dead.\nCheck the log file $log.\n\n");
  }
}

###########################################################################

sub Percent
#
# Print percentage
#
{
  my $p = 100 * $_[0];

  printf "\b\b\b\b%3.0f%%", $p;
}

sub Percent1
{
  my $p = 100 * $_[0];

  printf "\b\b\b\b\b\b%5.1f%%", $p;
}

sub Percent2
{
  my $p = 100 * $_[0];

  printf "\b\b\b\b\b\b%5.2f%%", $p;
}

#######################################################

sub Stamp
{
  local *CMD;
  open CMD, ">>command.log" or die;
  my @now = localtime();
  my $timestamp = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $now[5]+1900, $now[4]+1, $now[3], $now[2], $now[1], $now[0]);
  my @script = split /\//, $0;
  print CMD "$timestamp $script[-1] ", join(" ", @ARGV), "\n";
  close CMD;
}

sub sStamp
{
  my @now = localtime();
  my $timestamp = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $now[5]+1900, $now[4]+1, $now[3], $now[2], $now[1], $now[0]);
  my @script = split /\//, $0;
  return "$timestamp $script[-1] " . join(" ", @ARGV);
}

#######################################################

sub ReadASCII
#
# Reads columns from an ASCII file into an array
#
# @arr = ReadASCII($filename)
#
{
  my $file = shift;
  local *IN;
  open(IN, "$file") or die "Arrrgh: $!\n";
  my @in = ();
  while(<IN>)
  {
    my @arr = split;
    push @in, [@arr];
  }
  close(IN);
  return @in;
}


################################################

sub ReadFileList
#
# Read list of pathnames from an ASCII file
#
#   @list = ReadFileList($file_name);
#
# INPUT:
#
#   file_name - name of a file, containing list of pathnames
#
# OUTPUT:
#
#   list - list of pathnames
#
{
  my $file_name = $_[0];
  my @file_list = ();
  my $line;

  local *DATA;
  open(DATA, $file_name) 
    or die("\n*** Cannot open input file $file_name.\n\n");
  while($line = <DATA>) 
  {
    chomp $line;
    push(@file_list, $line) if $line !~ /^\s*$/;
  }
  close(DATA);
  return @file_list;
}

################################################


sub ReadLine
#
# Reads a line of text into an array, separator is space
#
#    @arr = ReadLine($line);
#
# INPUT:
#
#   line - a line of text
#
# OUTPUT
#
#   array - items from the line
#
{
  my $line = $_[0];
  my @arr = ();

  chop $line;
  $line =~ /( *)(.*)/;
  @arr = split(/ +/, $2);
  return @arr;
}


################################################################

sub ReadName
#
# Read a name from standard input
#
#   $answer = ReadName("$query", "$defaut");
#
# INPUT:
#
#   query - a question asked
#   default - default answer
#
# OUTPUT:
#
#   answer - answer typed; if only Enter hit, <default> is returned
#
# Queries and returns the answer. If return is typed, the default answer 
# is returned.
#
{
  my ($query, $default_answer) = @_;
  local *DIR;

  opendir(DIR, '.') or die $!;
  my @d = readdir(DIR);
  my $answer = Complete("$query [$default_answer]: ", \@d);
  if($answer eq "") {$answer = $default_answer;}
  if($answer eq "kill") {die "\nYou killed me, I'm dead, good bye.\n";}
  return $answer;
}

################################################

sub ReadYesNo
#
# Read a yes-or-no answer from standard input
#
#   $answer = ReadYesNo("$query", $defaut);
#
# INPUT:
#
#   query - a question asked
#   default - default answer (0/1)
#
# OUTPUT:
#
#   answer - 0 or 1
#
#
{
  my ($query, $default_answer) = @_;

  my $answer;
  my $def = ($default_answer) ? "y" : "n";
  do{
    print "$query (y/n) [$def]:";
    $answer = <STDIN>;
    chomp $answer;
  } until ($answer eq "" || $answer =~ /^[yn]$/);
  if($answer eq "") {$answer = $default_answer;}
#  if($answer eq "kill") {die "\nYou killed me, I'm dead, good bye.\n";}
  if($answer eq "y") {$answer = 1;}
  if($answer eq "n") {$answer = 0;}
  return $answer;
}

################################################################

sub Retrieve
#
# Retrieves parameters from the file
#
#  Retrieve($filename, $par1, $par2, ..., $parn);
#
# INPUT:
#
# $filname - name of the file (e.g. ".rc")
#
# OUTPUT
#
# $par<i> - paramter
#
{
  my $file_name = shift(@_);
  my ($par);

  if(-e $file_name)
  {
    local *RC;
    open(RC, $file_name) or die("Cannot open resource file $file_name.\n");
    my $i = 0;
    while($par = <RC>)
    {
      chop($par);
      $_[$i] = $par;
      $i++;
    }
  }
}

################################################################

sub Rmmin
#
# Removes minus signs breaking lines in QDP output from XSPEC.
# Lines are joined together.
#
{
  my $filename = shift;
  local (*IN, *OUT);
    
  my $tmpfile = "${filename}~";
  move($filename, $tmpfile) or die "$!\n";

  open(IN, "$tmpfile") or die "$!\n";
  open(OUT, ">$filename") or die "$!\n";
  while(my $line = <IN>)
  {
    chomp($line);
    if($line =~ /(.*)(-$)/)
      {print OUT "$1";}
    else
      {print OUT "$line\n";}
  }
  close(IN);
  close(OUT);
}

################################################################

sub Store
#
# Stores parameters in a file
#
#  Store($filename, $par1, $par2, ..., $parn);
#
# INPUT:
#
# $filname - name of the file (e.g. ".rc")
# $par<i> - paramter
#
{
  my $file_name = shift;

  local *RC;
  open(RC, ">$file_name") or die("Cannot open recource file $file_name.\n");
  foreach my $par (@_)
  {
    print RC "$par\n";
  }
  close(RC);
}

################################################################

sub WriteASCII
#
# Writes columns from an array into ASCII file
#
# WriteASCII(\@arr, $filename)
#
{
  my $data = $_[0];
  my $filename = $_[1];
  local *TMP;
  
  open(TMP, ">$filename") or die "$!\n";;
  foreach my $bin (@$data)
  {
    print TMP "@$bin\n";
  }
  close(TMP);
}

###################################################################

sub ReadFastaSequences
#
# my %S = ReadFastaSequences($file);
#
# Reads all sequences from a FASTA file into a hash 
#
#    seq_id => [$head, $sequence]
#
# The key 'seq_id' is the beginning of the comment line (without '>') until
# the first white space character. $head is the entire header line (without '>').
# These identifiers MUST be unique.
#
{
  my ($file, $verbouse) = @_;
  
  my %S = ();
  local *IN;
  open IN, $file or die "\nCannot open file $file\n";
  my $seq = '';
  my $head = '';
  my $id = '';
  while(my $line = <IN>)
  {
    Chomp($line);
    next if $line =~ /^;/;  # comment
    if($line =~ s/^>//)
    {
      if($verbouse)
      {
        my $mem = sprintf "%.3g", total_size(\%S) / 1e9;
        print "$line \(memory $mem GB)\n";
      }
      
      $S{$id} = [$head, $seq] if($id ne '');
      $head = $line;
      if($line =~ /^(\S+)\s+(\S.*)$/)
        {$id = $1}
      else
        {$id = $line}
      $seq = '';
    }
    else
      {$seq .= $line;}
  }
  $S{$id} = [$head, $seq] if($id ne '');
  return %S;
}

sub Chomp {$_[0] =~ s/(\n|\r\n|\r)$//}
###################################################################

sub WriteFastaSequences
{
  my ($file, $S, $comment) = @_;
  
  local *OUT;
  open OUT, ">$file" or die "\nCannot open output file $file\n";
  print OUT ";$comment\n" if defined $comment;
  for my $id (sort keys %$S)
  {
    my ($head, $seq) = @{$S->{$id}};
    print OUT ">$head\n";
    for my $i (0 .. int(length($seq) / 80))
      {print OUT substr($seq, $i * 80, 80), "\n";}
  }
  close OUT;
}


###############################################################

sub wc
{
  my $file = shift;
  my $count = `wc -l < $file`;
  die "wc failed: $?" if $?; 
  chomp $count ;
  return $count;
}

#########################################################

sub commify
#
# Format large integer number with commas, e.g. 10289371 => 10,289,371 
# 
{
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text
}

#########################################################

sub SmoothMovingMean
#
# Moving average with weighted window
#
{
  my ($p, $winsize, $outzeroes) = @_;

  $winsize = 1 unless defined $winsize;  
  my $N = nelem($p);
  
  my $lside = zeroes($winsize);
  my $rside = zeroes($winsize);
  unless($outzeroes)
  {
    $lside +=  at($p, 0);
    $rside += at($p, $N - 1);
  }
  my $P = append(append($lside, $p), $rside);

  # create weight window
  my $M = 2 * $winsize + 1;
  my $weight = zeroes($M);
  for my $i (1 .. $winsize + 1)
  {
    $weight($i - 1) .= $i;
    $weight($M - $i) .= $i;
  }
  my $norm = sum($weight);
  
  my $s = zeroes($p);
  for my $i (0 .. $N - 1)
  {
    my $j = $i + $winsize;
    $s($i) .= sum($P($j-$winsize:$j+$winsize) * $weight) / $norm;
  }
  return $s;
}

 

1;
