package Cluster;

our $VERSION = '0.1';

=head1 NAME 

Cluster -- Module to manage submissions to the cluster

=head1 SYNOPSIS


=head1 DESCRIPTION



=head1 AUTHOR

Chris Cole <christian@cole.name>

=head1 COPYRIGHT

Copyright 2007, Chris Cole.  All rights reserved.

This library is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut

use strict;
use warnings;
use Carp;

my $defaultSettings = 'source /grid.master/default/common/settings.sh &&';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(submit_cluster job_status);

#
## Submit jobs to the cluster
#
sub submit_cluster ($;$) {
   
   my ($cmd, $q) = @_;
   
   my $queue;
   if ($q eq '64') {
      $queue = '64bit.q';
   } elsif ($q eq '64pri') {
      $queue = '64bit-pri.q -l qname=64bit-pri.q'
   } else {
      $queue = '64bit.q';
   }
   my $qsub = "qsub -q $queue -cwd -S /usr/bin/perl";
   
   my $out = `$defaultSettings $qsub $cmd` or croak "ERROR - qsub failed for '$cmd'";
   my $id = (split /\s+/, $out)[2];
   if ($id !~ /^\d+$/) {
      croak "ERROR - qsub didn't return a job ID for '$cmd': $out";
   }
   return($id);
}

#
## Check the status of submitted jobs
#
sub job_status ($) {
   # Requires: a jobid as input
   # Returns:  0 if job not found
   #           1 if job is found and queuing
   #           2 if job is found and running
   #          -1 if job is in error state 
   
   my ($job) = @_ or return;
   
   chomp(my $me = `whoami`);
   
   my $cmd = "$defaultSettings qstat -u $me | grep $job";
   my $out = `$cmd`;
   
   return(0) if (!$out);
   
   my $state = (split /\s+/, $out)[4];
   
   if ($state =~ /E/) {
      # job has errored
      return(-1);
   } elsif ($state eq 'qw') {
      # job is queuing
      return(1);
   } else {
      # assume job is running
      return(2);
   }
}
