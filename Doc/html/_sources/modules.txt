.. modules:

.. _modules:

************
Perl modules
************

.. BootstrapDB_perlmod:

.. _BootstrapDB_perlmod:

**************
BootstrapDB.pm
**************

Set of routines to create/manipulate Sqlite databases for DE bootstrap.

.. Cluster_perlmod:

.. _Cluster_perlmod:

**********
Cluster.pm
**********

Module to manage submissions to the cluster

.. Distribution_perlmod:

.. _Distribution_perlmod:

***************
Distribution.pm
***************

Set of routines for distribution testing statistics.
		
.. GetOpt_perlmod:

.. _GetOpt_perlmod:

*********
GetOpt.pm
*********

GetOpt addon for the scripts

:param: `-psfile` - Output postscript file.

:param: `-cond` - Condition to use. In some scripts analysis can be restricted to one condition.

:param: `-rep` - Replicate number.

:param: `-nonzero` - Selects only data rows (genes) with at least *n* non-zero replicates.

:param: `-norm` - Expression data normalization. Input expression files contain raw counts. These are normalized on the fly, using one of the following methods:

  * none - no normalization
  * deseq - DESeq normalization (default)
  * tmm - trimmed mean of M values normalization
  * fpkm - approximate FPKM normalization (not identical to cuffdiff!)
  * totcount - total count
  * totcountm - total count with normalization factors scaled to the mean of 1

:param: `-type` - Type of raw data, as encoded in the counts file name (e.g. WT_raw.tsv is 'raw'). Could be 'lev' for leveled data, or some other, alternative prescription.

:param: `-exclude` - A list of replicates to exclude from analysis. The format of this string can be, e.g.
   
   * *1,3,5:* to exclude replicates 1, 3 and 5 from the first condition
   * *:1,2,3* to exclude replicates 1, 2 and 3 from the second condition
   * *1:3,4* to exclude replicates 1 from the first condition and 3 and 4 from the second condition  

:param: `-include` - A list of replicates to include to analysis. Format as above.

:param: `-clip` - Clip *n* top and n bottom replicates.

:param: `-genes` - Comma-delimited list of genes.

:param: `-genlist` - A file with a list of genes (first column).

:param: `-randrep, -randrep2` - Internal options.

:param: `-clean` - A switch indicating that only clean replicates should be used. Clean replicates are defined in :ref:`GRNASeq.pm <GRNASeq_perlmod>` module.

:param: `-multicor` - Multiple test correction to apply to tests/tools that return uncorrected raw p-values. The default value is *BH*.
 
:param: `-countsdir` - Path to the directory containing gene read count data. There should be two files, one for each conditions, with names C<cond>_raw.tsv, created by :ref:`combine_replicates.pl <combine_replicates_perldoc>`. If not specified, the value from *defs.dat* will be used.

:param: `-outlierdir` - A directory to store outlier data. If not specified, the value from *defs.dat* will be used.
 
:param: `-pileupdir` - A directory to store pileup data. If not specified, the value from *defs.dat* will be used.
 
:param: `-gff` - A GTF file to be used by many scripts. If not specified, the value from *defs.dat* will be used.
 
:param: `-clean` - A switch indicating that only clean replicates should be used. Clean replicates are defined in C<GRNASeq.pm> module.
 
:param: `-spclean` - A switch indicating that only spikein-clean replicates should be used. Spkein-clean replicates are defined in :ref:`GRNASeq.pm <GRNASeq_perlmod>` module.
 
.. GRNASeq_perlmod:

.. _GRNASeq_perlmod:

**********
GRNASeq.pm
**********


.. ReadExpressionFile_perldoc:

.. _ReadExpressionFile_perldoc:

------------------
ReadExpressionFile
------------------

Read data from expression file *$file*. Expression file should be 
tab-delimited with the first column containing gene/transcript name, and 
the remaining columns containing replicates. This script returns a 
two-dimensional piddle *$d* (genes in rows and replicates in columns) and 
one-dimensional piddles *$m* and *$s*, containing mean and standard deviation, 
respectively. The last value is a selection piddle::

	my ($d, $m, $s, $sel) = ReadExpressionFile($file, %options);

:param: `nonzero` - selects only rows (genes) with at least n non-zero replicates.
  
:param: `exclude [i1, i2, ...,in]` - excludes a list of replicates from data.

:param: `ranrep` - selects n replicates at random

.. ReadGenelist_perldoc:

.. _ReadGenelist_perldoc:

------------
ReadGenelist
------------

Reads a list of genes (first column) from expression file::

	my @g = ReadGeneList($file);

.. ReadTwoDETestFiles_perldoc:

.. _ReadTwoDETestFiles_perldoc:

------------------
ReadTwoDETestFiles
------------------

Reads two files with results of DE analysis. Each of these files should 
contain gene names in the first column, log2 fold changes in the second 
column and p-values in the third column.

Returned are: array of genes $genes, and two piddles with corresponding 
p-values, *$p1* and *$p2* and fold changes *$fc1*, *$fc2*:: 

	my ($genes, $p1, $p2, $fc1, $fc2) = ReadTwoDETestFiles($file1, $file2);

.. ReadDETestFile_perldoc:

.. _ReadDETestFile_perldoc:

--------------
ReadDETestFile
--------------

Reads a file with results of DE analysis. It should contain gene names in the 
first column, log2 fold changes in the second column and p-values in the third 
column. Returned is a hash: *gene => [fc, p]*::
	
	my %g = ReadDETestFile($file);


.. ReadSpikeAndDETestFiles_perldoc:

.. _ReadSpikeAndDETestFiles_perldoc:

-----------------------
ReadSpikeAndDETestFiles
-----------------------

Reads a DE file (fcp format) and a spike-in file (a list of good spikes with 
columns: spikein *@col1* and expected fold change *@col5*) and combines them 
together. Returns a list of gene (spike-in) names, spike-in significance 
(0 for not changing, 1 for changing) and p-values from the DE file::

	my ($genes, $sp, $p) = ReadSpikeAndDETestFiles($DEfile, $spikefile)


.. NormalizeData_perldoc:

.. _NormalizeData_perldoc:

-------------
NormalizeData
-------------

Normalizes data read using ReadExpressionFile: two-dimensional piddle with 
genes in rows and replicates in columns. Methods are 'none', 'deseq', 
'totalcount', 'tmm', 'fpkm'::

	my $nd = NormalizeData($d, $method);


.. GetNormalizedData_perldoc:

.. _GetNormalizedData_perldoc:

-----------------
GetNormalizedData
-----------------

Read data from expression file, select genes/replicates and normalize.

Returns the following:
  
	* $d - two-dimensional data piddle (rows = genes, columns = replicates)
	* $m, $s - mean and standard deviation of genes
	* $nrep - number of replicates (columns)
	* $ngen - number of genes
	* $genes - full list (array reference) of genes, before selection
	* $name - string to display in figures
	* $sel - selection of genes, can be used with $genes array
	* $f - normalizing factors

::

	my ($d, $m, $s, $nrep, $ngen, $genes, $name, $sel, $f) = GetNormalizedData(%options);

:param: `experiment` - not used yet

:param: `condition` - condition name (WT, Snf2...)

:param: `type` - type of data (raw, fake)

:param: `nonzero` - required number of non-zero replicates (-1 = all)

:param: `norm` - normalization (deseq, tmm, totalcount, none)

:param: `exclude` - exclude a list of replicates from data (comma separated)

:param: `clip` - clip n top and n bottom replicates  

.. PlotDistWithFits_perldoc:

.. _PlotDistWithFits_perldoc:

----------------
PlotDistWithFits
----------------

Fit data with theoretical distributions and plot histogram/KDE of data and 
best-fitting distributions::

	PlotDistWithFits($win, $pan, $d, $lab, %opt);

	* $win - window object (see PLGraphs)
	* $pan - panel number (1, 2, 3, ...)
	* $d - data (one-dimensional piddle)
	* $lab - label in the figure

:param: `charsize` - character size

:param: `hist` - plot histogram

:param: `kde` - plot KDE

:param: `onebin` - histogram in bins of 1

:param: `bestbin` - select best binning automatically

:param: `max` - maximum in y axis 

.. ReadTestfileInfo_perldoc:

.. _ReadTestfileInfo_perldoc:

----------------
ReadTestfileInfo
----------------

Reads DE tool information from a metafile *$file*::

	my %info = ReadTestfileInfo($file)

.. housekeeping_autodoc:

.. _housekeeping_autodoc:

***************
housekeeping.py
***************

.. automodule:: housekeeping

----------------------------------------
Custom Warnings and pretty times & dates
----------------------------------------

.. autofunction:: custom_formatwarning
.. autofunction:: timeAndDateStr
.. autofunction:: timeStr

--------------
Option Parsing
--------------

.. autofunction:: parse_required_options
.. autofunction:: parse_allowed_values
.. autofunction:: parse_option_type

-------
Logging
-------

.. autofunction:: write_log_header

-------------------------------
Temporary directories and files
-------------------------------

.. autofunction:: createNewTempdir

.. Powertests_perlmod:

.. _Powertests_perlmod:

*************
Powertests.pm
*************

.. rerr_perldoc:

.. _rerr_perldoc:

----
rerr
----
Calculates $x / ($x + $y) and its error::

	my ($r, $sr) = rerr($x, $y, $sx, $sy);

.. GetTotData_perldoc:

.. _GetTotData_perldoc:

----------
GetTotData
----------

Reads data from powertest file::
	
	my ($nrep, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetTotData($file)

.. GetFCData_perldoc:

.. _GetFCData_perldoc:

---------
GetFCData
---------

Reads data from powertest FC file::

  my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetTotData($file)

.. GetFCOneData_perldoc:

.. _GetFCOneData_perldoc:

------------
GetFCOneData
------------

Reads data from powertest FC file and combines positive and negative FC::

	my ($nrep, $fclim, $tp, $tn, $fp, $fn, $stp, $stn, $sfp, $sfn) = GetTotData($file)

.. CollectFC_perldoc:

.. _CollectFC_perldoc:

---------
CollectFC
---------

Collects unique fold changes from a piddle::

	my @F = CollectFC($fclim);


.. Tools_perlmod:

.. _Tools_perlmod:

********
Tools.pm
********

.. TreeView_perlmod:

.. _TreeView_perlmod:

***********
TreeView.pm
***********

A few subroutines supproting treeview and hierarchical clustering.

---------------------------
HierarchicalClusterProfiles
---------------------------
::

  HierarchicalClusterProfiles($profiles, $method, $dist);

------------
WriteCDTFile
------------
::

  WriteCDTFile($file, $tree, $names, $descriptions, $profiles);


------------
WriteGTRFile
------------
::

  WriteGTRFile($file, $tree);

-----------
CollectTree
-----------
::

  @nodes = CollectTree($tree, $node);

----
Dice
----
::

  @clustered = Dice($arrayref, \@nodes);

..	toctree::
	:maxdepth: 1

