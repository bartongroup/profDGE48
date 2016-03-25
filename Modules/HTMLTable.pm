package HTMLTable;
require Exporter;

use strict;
use warnings;

our $VERSION = 1.00;

our @ISA = ("Exporter");
our @EXPORT = qw(
  HTML_Header
  HTML_TableHead
  HTML_TableRow
  HTML_TableFinish
  HTML_BodyFinish
  $HTML_bkgcol
);

our $HTML_bkgcol = '#ffffff';

########################################################################

sub HTML_Header
{
  my ($fh, $head, $info) = @_;
  $head = '' if !defined $head;
  $info = '' if !defined $info;
    
  print $fh <<EOF;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<HTML>
<HEAD>
  <TITLE></TITLE>
  <LINK REL=STYLESHEET HREF="style.css" TYPE="text/css">
  <META http-equiv=Content-Type content="text/html; charset=utf-8">
  <META content="Nice little PERL script" name=GENERATOR>
 <style type="text/css">
  <!--
    td.head {
      font: 10pt arial, helvetica, sans-serif;
      color: #000000;
      padding: 0.1em 0em 0.1em 0.3em;
    }
    td.text {
      font: 9pt arial, helvetica, sans-serif;
      color: #000000;
      padding: 0.1em 0em 0.1em 0.3em;
    }
    td.seq {
      font: 9pt 'Lucida Console', Monaco, monospace;
      color: #000000;
      padding: 0.1em 0em 0.1em 0.3em;
    }
    H1 {
      font: 18pt arial, helvetica, sans-serif;
      font-weight: bold;
      color: #000000;
      padding: 0.1em 0em 0.1em 0.3em;
    }
    H2 {
      font: 14pt arial, helvetica, sans-serif;
      font-weight: bold;
      color: #000000;
      padding: 0.1em 0em 0.1em 0.3em;
    }
    P.text {
      font: 10pt arial, helvetica, sans-serif;
      color: #000000;
      padding: 0.1em 0em 0.1em 0.3em;
    }
    P.text2 {
      font: 12pt arial, helvetica, sans-serif;
      color: #000000;
      padding: 0.1em 0em 0.1em 0.3em;
    }
      
    a:link {color: #ff0030; text-decoration: none;}
    a:visited {color: #ff0030; text-decoration: none;}
    a:hover {color: #000000; text-decoration: underline;}
    a:active {color: #ff0000; text-decoration: underline;}
   -->
  </style>
</HEAD>

<BODY>

<div class="box">

<P style="font:24pt Verdana,Geneva,sans-serif;font-weight:bold">$head</P>

<P class="text2">$info</P>
EOF
}

########################################################################

sub HTML_TableHead
{
  my ($fh, $columns, $widths, $titles, $width, $bgcolor) = @_;
  $width = 1000 if !defined $width;
  $bgcolor ||= '#ffebcd';

  print $fh <<EOF;
<TABLE width="$width" cellSpacing="0" cellPadding="2" bgColor="#ffffff" border="0">
<TBODY>
<TR vAlign="top" align="center" bgColor="$bgcolor">
EOF

  my $i = 0;
  for my $col (@$columns)
  {
    my $w = (defined $widths && $widths->[$i] ne '') ? " width=\"$widths->[$i]\"" : '';
    my $t = (defined $titles && $titles->[$i] ne '') ? " title=\"$titles->[$i]\"" : '';
    print $fh "<TD $w class=\"head\" $t>$col</TD>\n";
    $i++;
  }
  
  print $fh "</TR>\n";
}

########################################################################

sub HTML_TableRow
{
  my ($fh, $columns, $widths, $class) = @_;

  $class ||= 'seq';
  
  print $fh "<TR vAlign=\"top\" bgColor=\"$HTML_bkgcol\">\n";
  my $i = 0;
  for my $col (@$columns)
  {
    my $w = (defined $widths && $widths->[$i] ne '') ? " width=$widths->[$i]" : '';
    print $fh "<TD $w class=\"$class\">$col</TD>\n";
    $i++;
  }
  
  print $fh "</TR>\n";
  $HTML_bkgcol = ($HTML_bkgcol eq '#ffffff') ? "#f0f0f0" : '#ffffff';  
}

########################################################################

sub HTML_TableFinish
{
  my $fh = shift;
  
  print $fh "</TBODY></TABLE>\n";
}

sub HTML_BodyFinish
{
  my $fh = shift;
  
  print $fh "</div></BODY>\n";
}





1;