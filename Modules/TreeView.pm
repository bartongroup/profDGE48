package TreeView;

=head1 NAME

CompBio::TreeView - a few subs supproting treeview and hierarchical clustering

=head1 SYNOPSIS

  HierarchicalClusterProfiles($profiles, $method, $dist);
  WriteCDTFile($file, $tree, $names, $descriptions, $profiles);
  WriteGTRFile($file, $tree);
  @nodes = CollectTree($tree, $node);
  @clustered = Dice($arrayref, \@nodes);

=head1 SUBROUTINES

=over 4

=item B<WriteCDTFile($file, $tree, $names, $descriptions, $profiles)>

Write CDT file for TreeView. $tree is created with Algorithm::Cluster::treecluster. Three last arguments are array references to names, descriptions and profiles (each profile is a piddle) in their original order before clustering. This subroutine will reorder these arrays according the $tree and print into $file in the format recognised by F<treeview>.

=item B<WriteGTRFile($file, $tree)>

Write GTR file for TreeView. $tree is created with Algorithm::Cluster::treecluster.

=back

=cut

require Exporter;
use strict;
use warnings;

use PDL;

our $VERSION = 0.1;

our @ISA = ("Exporter");
our @EXPORT = qw(
  WriteIDFile
  WriteCDTFile
  WriteGTRFile
  WriteNewickFile
  CollectTree
  Dice
  Newick
  HierarchicalClusterProfiles
);

###########################################################

sub WriteIDFile
{
  my ($file, $tree, $names) = @_;
  
  local *FILE;
  open FILE, ">$file" or die "\nCannot open output file $file\n";
  
  my @ids = CollectTree($tree);
  my @name = Dice($names, \@ids);
  
  for my $id (@name)
    {print FILE "$id\n";}
  close FILE;
}

###########################################################

sub WriteCDTFile
{
  my ($file, $tree, $names, $genes, $descriptions, $profiles) = @_;
  
  local *OUT;
  open OUT, ">$file" or die "\nCannot open CDF file $file for writing.\n";
  
  my @ids = CollectTree($tree);
  my @prof = Dice($profiles, \@ids);
  my @name = Dice($names, \@ids);
  my @desc = Dice($descriptions, \@ids);
  my @gene = Dice($genes, \@ids);

  my $proflen = nelem($profiles->[0]);
  print OUT "GID\tNAME\tGENE\tDESCRIPTION\tGWEIGHT";
  print OUT (map "\tP$_", (0 .. $proflen - 1)), "\n";
  
  print OUT "EWEIGHT\t\t\t\t1";
  print OUT "\t" x ($proflen - 1), "\n";
  
  for my $i (0 .. @ids - 1)
  {
    print OUT "PROT", $ids[$i], "\tgi:", $name[$i], "\t", $gene[$i], "\t", $desc[$i], "\t", 1;
    for my $n (0 .. $proflen - 1)
      {printf OUT "\t%6.3f", 2 * at($prof[$i], $n) - 1;}   # renormalize to -1, 1 for nice display in TreeView
    print OUT "\n";
  }
  close OUT;
}


###########################################################

sub WriteGTRFile
{
  my ($file, $tree, $rank) = @_;
  
  my %corr = rank_distances($tree) if $rank;
  
  local *OUT;
  open OUT, ">$file" or die "\nCannot open CDF file $file for writing.\n";

  print OUT "NODEID\tLEFT\tRIGHT\tCORRELATION\n";

  my $N = $tree->length;  
  for my $i (0 .. $N - 1)
  {
    my $node = $tree->get($i);
    my $D = ($rank) ? $corr{$node->distance} : 1 - $node->distance; 
    printf OUT "%s\t%s\t%s\t%9.5f\n",nnum(-1-$i),nnum($node->left),nnum($node->right), $D;
 }
}

sub nnum
{
  my $n = shift;
  ($n >= 0) ? sprintf "PROT%d", $n : sprintf "NODE%d", -$n;
}

###########################################################

sub WriteNewickFile
{
  my ($file, $tree, $ids) = @_;
  
  local *OUT;
  open OUT, ">$file" or die "\nCannot open Newick file $file for writing.\n";
  my ($S) = Newick($tree, undef, $ids);
  print OUT "$S\n";
  close OUT
}

###########################################################

sub CollectTree
{
# returns all nodes of the subtree of the given node

  my ($tree, $n) = @_;
  
  $n = -$tree->length if !defined $n;  # top node
  return ($n) if $n >= 0;
  my $node = $tree->get(-1 - $n);
  my @left_list = CollectTree($tree, $node->left);
  my @right_list = CollectTree($tree, $node->right);
  return (@left_list, @right_list);
}

###########################################################

sub Newick
#
# Returns string representing tree in Newick format
#
{
  my ($tree, $n, $names) = @_;
  
  $n = -$tree->length if !defined $n;  # top node
  if($n >= 0)
  {
    my $leaf = (defined $names) ? $names->[$n] : $n;
    $leaf =~ s/\"//g;
    $leaf =~ s/[\s\,\:\;]/_/g;
    return ($leaf, 0);
  }
  my $node = $tree->get(-1 - $n);
  my $distLR = $node->distance;        # distance from node to the leaves
  my ($L, $distL) = Newick($tree, $node->left, $names);
  my ($R, $distR) = Newick($tree, $node->right, $names);
  my $lenL = sprintf "%.3g", $distLR - $distL;  # lenght of left branch
  my $lenR = sprintf "%.3g", $distLR - $distR;  # length of right branch
  my $str = "($L:$lenL,$R:$lenR)";
  
  return ($str, $distLR);
}

###########################################################

sub Dice
{
  # Reorder array according to index, similar to PDL's dice
  my ($arrayref, $index) = @_;
  
  my @s = ();
  for my $i (@$index)
    {push @s, $arrayref->[$i];}
  return @s;
}

###########################################################

sub HierarchicalClusterProfiles
{
  my ($profiles, $method, $dist, $scale) = @_;
  
  $method = 'a' if !defined $method;
  $dist = 'c' if !defined $dist;
  
  my @data = ();
  my @mask = ();
  for my $p (@{$profiles})
  {
    my @prof = list $p;   # convert PDL into Perl list
    push @data, [@prof];
    push @mask, [(1) x nelem($p)];
  }
  my @weight = @{$mask[0]};
  
  my %params = (
    transpose => 0,
    method => $method,
    dist => $dist,
    data => \@data,
    mask => \@mask,
    weight => \@weight,
  );

  my $tree = Algorithm::Cluster::treecluster(%params);
  $tree->scale if $scale;  # does it do anything?
  return $tree;
}

###########################################################

sub rank_distances
{
  # rank and scale correlations for nice display in TreeView

  my $tree = shift;
  my %corr = ();
  my $N = $tree->length;  
  for my $i (0 .. $N - 1)
  {
    my $node = $tree->get($i);
    $corr{$node->distance} = 1;
  }
  my $num = scalar keys %corr;
  my $i = 0;
  for my $c (reverse sort {$a <=> $b} keys %corr)   # establish ranking in correlations
  {
    my $rankcorr = $i / ($num - 1);
    $corr{$c} = $rankcorr;
    $i++;
  }
  return %corr;
}






