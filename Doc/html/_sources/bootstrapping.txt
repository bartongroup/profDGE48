.. bootstrapping:

.. _bootstrapping:

***********************
Bootstrapping DGE calls
***********************

.. add_gene_name_column_perldoc:
		
.. _add_gene_name_column_perldoc:

=======================
add_gene_name_column.pl
=======================

add_gene_name_column.pl - given a tab delimited file with Ensembl IDs in the 
first col, add a gene name::

	add_gene_name_column.pl --in <file> --species <name> [--feature-type <string>] [--desc|--no-desc] [--coords|--no-coords] [--out <file>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help] [--version]

:param: `--in` - Input tab-delimited file. 1st column must be an ensembl ID.

:param: `--species` - Species name.

:param: `--feature-type` - Specify whether the Ensembl IDs correspond to genes or transcripts. [default: Genes]

:param: `--coords|--no-coords` - Toggle whether to add gene description as well as gene name. [default: off]

:param: `--out` - Output filename. [default: ensembl_annotated.csv]

:param: `--verbose|--no-verbose` - Toggle verbosity. [default: on]

:param: `--debug|--no-debug` - Toggle debugging output. [default: off]

:param: `--version` - Only print out the version information before exiting.

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

Given a tab-delimited list with Ensembl gene IDs in the first column, add gene 
names (and optionally descriptions and coordinates) to the list. The additional 
information is placed as the last set of columns in the data, leaving the 
original data at the same columns. If the file has a header line, the first 
column must be called 'GeneID' otherwise the formatting will be wrong.

.. combine_bs_pvalues_perldoc:

.. _combine_bs_pvalues_perldoc:

=====================
combine_bs_pvalues.pl
=====================

Calculates median of p-values across bootstrap runs. The input is all sqlite 
.db files created by DE bootstrap (for a range of replicate numbers). 
The script will create a single .tsv file per replicate number, containing 
median p-values. This script is a wrapper around :ref:`one_bs_pvalue.pl <one_bs_pvalue_perldoc>`, and
it runs it on the cluster::

	combine_bs_pvalues.pl -dbdir=./bs_results -dbform=de_ttest_rep%d.db

:param: `-dbdir` - Path to the directory containing .db files.

:param: `-dbform` - Format of .db file name, e.g. de_test_rep%d.db, for files de_test_rep2.db, de_test_rep3.db, ..., de_test_rep30.db. The replicates will be cycled between 2 and 48, non-existing files ignored.

:param: `-script` - This should be the full path to the :ref:`one_bs_pvalue.pl <one_bs_pvalue_perldoc>` script that this wrapper runs.

:param: `-logdir` - Path to the directory where all intermediate files are stored. (Default: ./de_tests_logs.)

:param: `-outform` - Format for output .tsv files.

:param: `-queue` - Cluster queue name or names (comma-separated).

:param: `-help` - Brief help.

:param: `--man` - Full manpage of program.

.. generic_wrapper_autodoc:

.. _generic_wrapper_autodoc:

******************
generic_wrapper.py
******************

.. automodule:: generic_wrapper

=============================
generic_wrapper.py Functions:
=============================

----------------------------
Experiment structure parsing
----------------------------

.. autofunction:: parse_experiment_data
.. autofunction:: runStructure

-------------
Bootstrapping
-------------

.. autofunction:: clusterBootstrap
.. autofunction:: clusterDistribute
.. autofunction:: rep_subselect
.. autofunction:: runClusterBootstrap

-------------
Data Routines
-------------

.. autofunction:: get_gene_counts
.. autofunction:: read_expression_data
.. autofunction:: run_agnc
.. autofunction:: run_gbg_cluster
.. autofunction:: run_gbg

-------------
Normalization
-------------

.. autofunction:: deseq_normalization
.. autofunction:: normalize_data
.. autofunction:: parseFlagstatOut
.. autofunction:: run_flagstat
.. autofunction:: run_flagstat_cluster

.. genwrap_Rinterface:

.. _genwrap_Rinterface:

--------------
Interface to R
--------------

.. autofunction:: write_R_data_file
.. autofunction:: runRDE

--------------------
General housekeeping
--------------------

.. autofunction:: concatSqliteFiles
.. autofunction:: monitorDrmaaJobs
.. autofunction:: results2sqlite
.. autofunction:: writeCheckpoint

.. grid_launcher_DE_perldoc:

.. _grid_launcher_DE_perldoc:

===================
grid_launcher_DE.pl
===================

This is a wrapper script that divides gene counts into small chunks (default 
30 genes each) and runs DE scripts on the cluster, using DRMAA. The results are 
then aggregated into one file::

  grid_launcher_DE.pl -script=./one_bootstrap_ET93.pl -logdir=./boot_output -outfile=boot_result.tsv -options nboot=100000
  grid_launcher_DE.pl -script=./one_permutest.pl -logdir=./perm_output -outfile=perm_result.tsv
  
:param: `-script` - File name (with absolute path) containing the script to run DE. The script should have two mandatory options, *<-infile>* and *<-outfile>* and any other additional options. Input file is a chunk of gene counts in individual replicates and both conditions. There are three lines per each gene: gene name, counts for condition one (tab-delimited) and counts for conditions two. Such file is created by this script and passed on to the DE script. Output file countains resulting DE calles. There should be at least five tab-delimited columns in it: gene name, log2 fold change, p-value, mean of condition 1, mean of condition 2. Individual output files are aggregated into the final output file. Any additional parameters can be included. These are passed on by *<-options>* option. DE scripts available at this moment are :ref:`one_bootstrap_ET93.pl <one_bootstrap_ET93_perldoc>` and :ref:`one_permutest.pl <one_permutest_perldoc>`.

:param: `-logdir` - Path to the directory where all intermediate files are stored.

:param: `-outfile` - Aggregated output file containing DE results.

:param: `-batchsize` - Number of genes in one batch. (Default: 30).

:param: `-options` - All other options required by the DE script, for example, number if iterations. See synopsis for an example. 

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.


.. grid_launcher_powertest_cuffdiff_perldoc:

.. _grid_launcher_powertest_cuffdiff_perldoc:
		
===================================
grid_launcher_powertest_cuffdiff.pl
===================================

This script launches cuffdiff power tests on the cluster. Cuffdiff is run 
*<nsim>* times for increasing number of replicates, starting from 2 and 
ending at *<maxn>*. Fold changes and p-value is calculated for each gene, for 
the given number of replicates across all simulations. These p-values can be 
later used to compute power of the test as a function of number of replicates::

	grid_launcher_powertest_cuffdiff.pl -datapath $g/mapping/genome_biol_reps -annotation Saccharomyces_cerevisiae.EF4.68.gtf -genome Scerevisiae68_ERCC92.fasta -genlist genlist.tsv -nsim=30 -maxn=40 -cuffdiff_script=./run_cuffdiff.pl
    
:param: `-cuffdiff_script` - File name (with absolute path) containing the script :ref:`run_cuffdiff.pl <run_cuffdiff_perldoc>`.

:param: `-datapath` - Path to the data directory of the experiment. Each sub-directory in this will be treated as a separate condition in the experiment (with the name of the condition matching the name of the directory), and each .bam file in each directory is a replicate for that condition.

:param: `-annotation` - Path to the .gff feature annotation file for the data. This file should match the feature you are counting the RNA-Seq expression for.

:param: `-genome` - Path to the FASTA file with reference genome. This is used by cuffdiff to run its bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates. You should have write access to the directory with this file as cufflinks will attempt to write an index file there. It might be a good idea to copy genome file (or a link) into your local directory.

:param: `-genlist` - Path to the file containing a list of all genes of interest. It can contain multiple columns, the first column should contain (case-insensitive) gene names.

:param: `-maxn` - Maximum number of replicates to consider. The scrip will iterate the DE test from 2 replicates up to this number. It should be smaller than the total number of replicates available!

:param: `-nsim` - Number of bootstrap runs for the given number of replicates. Median p-values will be calculated over these runs. The default value is 30.

:param: `-cuffdir` - Path to a directory where all cuffdiff output files will be stored. It will be divided into subdirectories for each number of replicates and each iteration. These are all intermediate files, so they can be deleted when you are satisfied with the result.

:param: `-randcond` - If specified, cuffdiff will be performed on randomly selected replicates from the same condition, defined here.

:param: `-queue` - Cluster queue name or names (comma-delimited).

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.
		
.. grid_launcher_powertest_perldoc:

.. _grid_launcher_powertest_perldoc:

==========================
grid_launcher_powertest.pl
==========================
		
Grid launcher for power tests. This script launches differential expression 
power tests on the cluster. The DE test of choice (see :ref:`simple_de_test.pl <simple_de_test_perldoc>` 
for details) is run *<nsim>* times for increasing number of replicates, 
starting from 2 and ending at *<maxn>*. P-values and fold changes are stored 
in sqlite databases, to be later processed by :ref:`combine_bs_pvalues.pl <combine_bs_pvalues_perldoc>` 
and :ref:`make_powerstats_db.pl <make_powerstats_db_perldoc>`. Note: this 
script is for simple tests (t-test, log t-test, MW). It does the same job 
as :ref:`generic_wrapper.py <generic_wrapper_autodoc>` does for tools like 
*edgeR* or *DEseq*. Requires script :ref:`one_powertest_simple.pl <one_powertest_simple_perldoc>`::

	grid_launcher_powertest.pl -test=st
	grid_launcher_powertest.pl -test=lt -nsim=1000 -logdir=./lt_test -options norm=deseq
  
:param: `-script` - File name (with absolute path) containing the script :ref:`one_powertest_simple.pl <one_powertest_simple_perldoc>`. This is a wrapper around the actual DE test script C<simple_de_test.pl>. If not specified, the dafult location, *<./scripts/one_powertest_simple.pl>*, will be used.

:param: `-logdir` - Path to the directory where all STDOUT and STDERR files from individual cluster jobs are stored. (Default: *<de_tests_logs>*).

:param: `-dbdir` - Path to the directory where all results are saved. These are Sqlite databases, one file per number of replicates. (Default: *<de_tests_db>*).

:param: `-test=I<[t|lt|st|mw|ks]>` - Type of differential expression test.

	* t - standard t-test for two samples with identical variances
	* lt - logarithmic t-test using approximately normal statistic Z = (log m2 - log m1) / sqrt(((s1/m1)**2)/n1 + ((s2/m2)**2)/n2), where m1, m2 are sample means and s1, s2 are sample standard deviations
	* st - shrinkage variance test of `Opgen-Rhein and Strimmer (2007) <http://www.degruyter.com/view/j/sagmb.2007.6.1/sagmb.2007.6.1.1252/sagmb.2007.6.1.1252.xml>`_
	* mw - Mann-Whitney test
	* ks - Kolmogorov-Smirnov test; **warning** - KS test does not work properly for discrete data!
  
:param: `-maxn` - Maximum number of replicates to consider. The scrip will iterate the DE test from 2 replicates up to this number. It should be smaller than the total number of replicates available!

:param: `-nsim` - Number of bootstrap runs for the given number of replicates. Median p-values will be calculated over these runs. (Default: 30).

:param: `-spike` - If specified, spike-ins will added to existing genes and spike-in-clean replicates will be selected.

:param: `-norm` - Normalization to be used (*deseq*, *totcountm*, *spikein*). (Default: *deseq*).

:param: -options` - Any additional options passed to the executed script (see Synopsis for usage).

:param: -queue` - Cluster queue name(s), comma-delimited.

:param: -randcond` - If specified, test will be performed on randomly selected replicates from the same condition, defined here.

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.		
		
.. group_by_gene_perldoc:

.. _group_by_gene_perldoc:

================
group_by_gene.pl
================

Essentially, this script is a wrapper for `htseq-count <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_, 
but makes it useable with BAM files. `htseq-count <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_ 
only works with SAM files sorted by readID 
whereas most other utilities require BAM files sorted by chromosomal position. 
This script takes a BAM file converts it into an appropriately sorted SAM file 
and runs `htseq-count <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_ 
on it. As SAM files can be pretty large and are only 
required temporarily here, the SAM file is stored on the TMPDIR path and 
deleted. Thereby, saving precious GPFS disk space. This behaviour can be c
ontrolled with the *<--local|--no-local>* switch. 

`htseq-count <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_ and 
`samtools <http://www.htslib.org/>`_ are required and must be in the default paths or 
provided with the *--samtools* and *--htseq* arguments.The script is optimised for 
the cluster so use the *<ram>* and *<local_free> I<SGE>* options to ensure the 
jobs don't run out of local disk space nor RAM::

	qsub -q 64bit-pri.q -l ram=4G -l local_free=70G group_by_gene.pl --bam my.bam --gtf my.gtf --out my.out

:param: `--bam` - Input BAM file.

:param: `--gtf` - Input GTF file.

:param: `--tool` - Select which tool to use for sorting (samtools or picard). [default: samtools]

:param: `--method` - Specify which method to use for aggregating reads to genes. [default: union]

:param: `--feature` - Specify which feature to aggregate read on. [default: gene_id]

:param: `--out` - Output filename. [default: gene.out]

:param: `--samtools` - Give path to samtools executable [default: /sw/samtools-0.1.18/samtools]

:param: `--htseq` - Give path to samtools executable [default: /sw/bin/htseq-count]

:param: `--local|--no-local` - Toggle whether to use /local for tmp files, otherwise use current dir. [default --local]

:param: `--verbose|--no-verbose` - Toggle verbosity. [default:none]

:param: `--debug|--no-debug` - Toggle debugging output. [default:none]

:param: `--version` - Print version info and exit

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. make_cuffdiff_db_perldoc:

.. _make_cuffdiff_db_perldoc:

===================
make_cuffdiff_db.pl
===================

Creates databases from cuffdiff bootstrap results. You need to run 
:ref:`grid_launcher_powertest_cuffdiff.pl <grid_launcher_powertest_cuffdiff_perldoc>`
first::

	make_cuffdiff_db.pl -cuffdir=cuffdiff_powertests -dbdir=de_tests_db -genlist=genes.lst   

:param: `-cuffdir` - Path to the directory where cuffdiff bootstrap results are stored.

:param: `-dbdir` - Path to the directory where all results are saved. These are Sqlite databases, one file per number of replicates.

:param: `-genlist` - File with gene names in the first column.

:param: `-sig>=[p|q]` - Whether to use p-value or q-value from the cuffdiff's output. Default is C<p>.

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.
		
.. one_bootstrap_ET93_perldoc:

.. _one_bootstrap_ET93_perldoc:

=====================
one_bootstrap_ET93.pl
=====================

Bootstrap test for two samples, to be run by grid launcher. Following 
`Efron & Tibshirani (1993), "An introduction to the bootstrap", p.220-224 
<http://www.hms.harvard.edu/bss/neuro/bornlab/nb204/statistics/bootstrap.pdf>`_.
To be used by :ref:`grid_launcher_DE.pl <grid_launcher_DE_perldoc>` and it is not 
supposed to be ran on its own.

:param: `-infile` - Name of the input file.

:param: `-outfile` - Name of the output file.

:param: `-nboot` - Number of bootstraps. The default value is 5000000.

.. one_bs_pvalue_perldoc:

.. _one_bs_pvalue_perldoc:

================
one_bs_pvalue.pl
================

Combines bootstrap results. To be used by 
:ref:`combine_bs_pvalues.pl <combine_bs_pvalues_perldoc>`. This script is not 
supposed to be run directly.

:param: `-dbfile` - Database file to process.

:param: `-outfile` - Output file.

:param: `-logit` - If this option is specified, log2  will be calculated of all fold-changes. However, we expect .db files to contain log2-fold-changes, so this is probably not necessary.

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.

.. one_permutest_perldoc:

.. _one_permutest_perldoc:

================
one_permutest.pl
================

Permutation test for two samples. To be used by :ref:`grid_launcher_DE.pl <grid_launcher_DE_perldoc>`.

:param: `-infile` - Input file.

:param: `-outfile` - Output file.

:param: `-nboot` - Number of bootstraps. The default number is 5000000.

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.

.. one_powertest_simple_perldoc:

.. _one_powertest_simple_perldoc:

=======================
one_powertest_simple.pl
=======================

Script to perform DE power test, used by 
:ref:`grid_launcher_powertest.pl <grid_launcher_powertest_perldoc>`. This is 
essentially a wrapper around :ref:`simple_de_test.pl <simple_de_test_perldoc>`. 
This script is not supposed to be ran on its own.

:param: `-test [t|lt|st|mw|ks]` - Type of differential expression test.

	* t - standard t-test for two samples with identical variances
	* lt - log-ratio t-test using approximately normal statistic Z = (log m2 - log m1) / sqrt(((s1/m1)**2)/n1 + ((s2/m2)**2)/n2), where m1, m2 are sample means and s1, s2 are sample standard deviations
	* st - shrinkage variance test of Opgen-Rhein and Strimmer (2007)
	* mw - Mann-Whitney test
	* ks - Kolmogorov-Smirnov test; B<warning> - KS test does not work properly for discrete data!
  
:param: `-nrep` - Number of replicates to run the test on.
  
:param: `-randcond` - If specified, test will be performed on randomly selected replicates from the same condition, defined here.

:param: `-nsim` - Number of bootstraps.
 
:param: `-dbfile` - Name of the sqlite .db file to create. 

:param: `-norm` - Expression data normalization. Input expression files contain raw counts. These are normalized on the fly, using one of the following methods:

	* none - no normalization
	* deseq - DESeq normalization (default)
	* tmm - trimmed mean of M values normalization
	* fpkm - approximate FPKM normalization (not identical to cuffdiff!)
	* totcount - total count
	* totcountm - total count with normalization factors scaled to the mean of 1

.. qsub_bootstraps_withexcludelist_shdoc:

.. _qsub_bootstraps_withexcludelist_shdoc:

==================================
qsub_bootstraps_withexcludelist.sh
==================================

Wrapper for running multiple instances of :ref:`generic_wrapper.py <generic_wrapper_autodoc>`
with different numbers of bootstrap iterations in each instance, and includes 
the use of an exclusion list for excluding specific replicates from the calculations.

.. qsub_bootstraps_withoutexludelist_shdoc:

.. _qsub_bootstraps_withoutexludelist_shdoc:

====================================
qsub_bootstraps_withoutexludelist.sh
====================================

Wrapper for running multiple instances of :ref:`generic_wrapper.py <generic_wrapper_autodoc>`
with different numbers of bootstrap iterations in each instance using the full 
set of replicates.
		
.. qsub_WTfdr_bootstraps_shdoc:

.. _qsub_WTfdr_bootstraps_shdoc:

========================
qsub_WTfdr_bootstraps.sh
========================

Wrapper for running multiple instances of :ref:`generic_wrapper.py <generic_wrapper_autodoc>`
with different numbers of bootstrap iterations in each instance using just one
condition of data. This was used to explore how well tool control their False
Discovery Rate.
	
.. simple_de_test_perldoc:

.. _simple_de_test_perldoc:
		
=================
simple_de_test.pl
=================

Simple differential expression test. Gene expression data should be stored in 
two tables, one for each condition. These can be created by the script 
:ref:`combine_replicates.pl <combine_replicates_perldoc>`. This script can 
perform a few flavours of t-test, Mann-Whitney test and Kolmogorov-Smirnov 
test. See option *-test* for details. Results are stored in a tab-delimited 
output file. Columns are:

	* gene name
	* log2 fold change
	* p-value
	* mean for first condition
	* mean for second condition

  simple_de_test.pl -norm=deseq -clean -test=mw

:param: `-test>=I<[t|lt|st|mw|ks]>` - Type of differential expression test.

	* t - standard t-test for two samples with identical variances
	* lt - logarithmic t-test using approximately normal statistic Z = (log m2 - log m1) / sqrt(((s1/m1)**2)/n1 + ((s2/m2)**2)/n2), where m1, m2 are sample means and s1, s2 are sample standard deviations
	* st - shrinkage variance test of `Opgen-Rhein and Strimmer (2007) <http://www.degruyter.com/view/j/sagmb.2007.6.1/sagmb.2007.6.1.1252/sagmb.2007.6.1.1252.xml>`_
	* mw - Mann-Whitney test
	* ks - Kolmogorov-Smirnov test; **warning** - KS test does not work properly for discrete data!
  
:param: `-outfile` - Output file. If not defined, a default name will be build using condition names and type of test.

:param: `-norm` - Expression data normalization. Input expression files contain raw counts. These are normalized on the fly, using one of the following methods:

	* none - no normalization
	* deseq - DESeq normalization (default)
	* tmm - trimmed mean of M values normalization
	* fpkm - approximate FPKM normalization (not identical to cuffdiff!)
	* totcount - total count
	* totcountm - total count with normalization factors scaled to the mean of 1

:param: `-clean` - A switch indicating that only clean replicates should be used. Clean replicates are defined in :ref:`GRNASeq.pm <GRNASeq_perlmod>` module.

:param: `-defsfile` - Path to defs.dat configuration file.

:param: `-randcond` - If specified, test will be performed on randomly selected replicates from the same condition, defined here.

..	toctree::
	:maxdepth: 1
