.. badreps:

.. _badreps:

****************************
Bad replicate identification
****************************

Scripts used to identify "bad" replicates in the experiemnt.

.. bad_replicates_perldoc:

.. _bad_replicates_perldoc:

=================
bad_replicates.pl
=================

Creates a combined figure to find bad replicates. Adds bad replicate numbers, as defined in C<defs.dat> file. Requires pileup statistics created with C<pileup_var.pl>::

	bad_replicates.pl

:param: `-limit` - Significance limit for outlier selection. This is a limit of number of standard deviations from the trimmed mean. (default: 5)

:param: `-ntrim` - Number of data points to remove on each side for the trimmed mean. (default: 3).

:param: `-off` - Offset for plotting bad replicate numbers above/below data points.

:param: `-psfile` - Name of the postscript file to redirect output to.

.. compare_replicates_perldoc:

.. _compare_replicates_perldoc:

=====================
compare_replicates.pl
=====================

Show various comparisons between replicates in a given condition::
	
	compare_replicates.pl -graph=all

There are several other options in this script, used for exploratory analysis
.
:param: `-graph [all|2]` - Type of plot to create.

	* *all* - heat plots with all replicates
	* *2* - two replicates (or selection of pairs)

:param: `-rep1` - First replicate to compare.

:param: `-rep2` - Second replicate to compare.

:param: `-cond` - Condition name.

:param: `-cond [WT|Snf2]` - Condition name.

:param: `-psfile` - Output postscript file.

..	toctree::
	:maxdepth: 1
