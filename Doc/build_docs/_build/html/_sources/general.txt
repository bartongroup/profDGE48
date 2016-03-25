.. General:

.. _General:

*************************
General scripts and files
*************************

General scripts and configuration files essential for the analysis.

.. combine_replicates_perldoc:

.. _combine_replicates_perldoc:

*********************
combine_replicates.pl
*********************

Combine counts from individual condition/replicate/lane files into multicolumn 
files. Genes in all files **must** be in the same order. Warning: replicate and 
lane number and file name templates hardcoded in the script.::

  combine_replicates.pl -indir=<input dir> -outdir=<output dir> -bioreps
  
:param: `-bioreps` -  tells the script to do biological replicates instead of lanes.

.. de_tests_configfile:

.. _de_tests_configfile:

************
de_tests.txt
************

Configuration file describing all tests and tools for the powertests. Entries 
should include:

  * name - name of the test to be displayed
  * repflie - format for bootstrap results (an FCP file) for all replicates (%d stands for replicate number)
  * sp_repfile = for spike-in bootstrap results (FCP file)
  * fullfile = result (FCP file) for full clean replicate set
  * adjustment = multiple test adjustment used by the tool to obtain p-values (none, BH or HB)
  * FCsig = optional sign correction for fold-change (+1 or -1)
  
 FCP files are tab-delimited with (at least) three columns: geneid, log2 
 foldchange and p-value.

.. de_tests_samecond_configfile:

.. _de_tests_samecond_configfile:

*********************
de_tests_samecond.txt
*********************

Configuration file describing all tests and tools for the null (same condition)
FDR tests. Entries should include:

  * name - name of the test to be displayed
  * repflie - format for bootstrap results (an FCP file) for all replicates (%d stands for replicate number)
  * sp_repfile = for spike-in bootstrap results (FCP file)
  * fullfile = result (FCP file) for full clean replicate set
  * adjustment = multiple test adjustment used by the tool to obtain p-values (none, BH or HB)
  * FCsig = optional sign correction for fold-change (+1 or -1)
  
 FCP files are tab-delimited with (at least) three columns: geneid, log2 
 foldchange and p-value.

 .. defs_dat_configfile:

.. _defs_dat_configfile:

********
defs.dat
********

Simple text configuration file describing file locations, clean replicates and 
condition IDs, for example.::

	topdir = rootpath
	countsdir = combined_counts
	outlierdir = outliers
	pileupdir = pileup
	bootstrapdir = bootstrapscat 
	gfffile = Saccharomyces_cerevisiae.EF4.68.gtf
	clean = 21,22,25,28,34,36:6,13,25,35
	spclean = 21,22,25,28,31,34,36:3,6,10,12,13,14,15,17,21,23,24,25,27,29,30,31,33,35,36,37,41,44
	conditions = WT,Snf2

.. exclude_badreps_lst_configfile:

.. _exclude_badreps_lst_configfile:

*******************
exclude_badreps.lst
*******************

List of replicate filenames (one per line) to be excluded from bootstrapping 
runs by :ref:`generic_wrapper.py <generic_wrapper_autodoc>`.

.. genlist_configfile:

.. _genlist_configfile:

***********
genlist.tsv
***********

Simple text file listing (one per line) the SGD gene names for Saccharomyces
cerevisiae.
