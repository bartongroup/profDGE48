.. plotting:

.. _plotting:

*****************
Plotting routines
*****************

.. cluster_replicates_pv_Rdoc:

.. _cluster_replicates_pv_Rdoc:

***********************
cluster_replicates_pv.R
***********************

This short Rscript uses the R package `pvclust <http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/pvclust/>`
to hierarchically cluster the experiment replicates using correlation distance
and complete linkage. The package boostraps the clustering process to assign 
significances for the structure of the observed clustering. This script uses
1000 bootstraps and marks cluster that are not observed in at least 95% of the 
boostraps.

.. cluster_tools_pv_Rdoc:

.. _cluster_tools_pv_Rdoc:

******************
cluster_tools_pv.R
******************

This short Rscript uses the R package `pvclust <http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/pvclust/>`
to hierarchically cluster the DGE tools tested here using correlation distance
and complete linkage. The package boostraps the clustering process to assign 
significances for the structure of the observed clustering. This script uses
1000 bootstraps and marks cluster that are not observed in at least 95% of the 
boostraps.

.. compare_null_perldoc:

.. _compare_null_perldoc:

***************
compare_null.pl
***************

Compare results from null test - FPR bootstraps for the same condition. 
Requires a test description file with all details of the bootstraps 
(*de_tests_samecond.txt* is default).

  compare_null.ps -psfile=null.ps

:param: `-testinfofile` - A file with DE tools metadata (for the 'same condition' bootstraps). See :ref:`plot_fcp.pl <plot_fcp_perldoc>` for details. Default is *de_tests_samecond.txt*.

:param: `-dir` - Directory where the results from :ref:`make_powerstats_db.pl <make_powerstats_db_perldoc>` are stored. The default value is *powerstats_same_db*.

:param: `-tests` - A comma-delimited list of tests/tools to be included in the plot. The default value is '*lt,bayseq,cuffdiff,degseq,deseq1,deseq2,ebseq,edger,limma,noiseq,samseq,poissonseq*'.

:param: `-graph` - Type of graph to plot. The default value is 'all'. Do not change!

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. compare_powerstats_db_perldoc:

.. _compare_powerstats_db_perldoc:

************************
compare_powerstats_db.pl
************************

Script to create figure 2.

  compare_powerstats_db.pl     
    
:param: `-testinfofile` - A file with DE tools metadata. See :ref:`plot_fcp.pl <plot_fcp_perldoc>` for details. Default is *de_tests.txt*.

:param: `-dir` - Directory where the results from :ref:`make_powerstats_db.pl <make_powerstats_db_perldoc>` are stored. The default value is *powerstats_db_ref*.

:param: `-truetest` - A test to be used as a 'gold standard'. The default value (as used in the Paper) is '*self*', that is the full clean replicate test for the same tool.

:param: `-tests` - A comma-delimited list of tests/tools to be included in the plot. The default value is '*lt,bayseq,cuffdiff,degseq,deseq1,deseq2,ebseq,edger,limma,noiseq,samseq,poissonseq*'.

:param: `-minnerp` - Minumum number of replicates to include in the plots. Default is 3.

:param: `-nrep` - Number of replicates to create the plot.

:param: `-graph` - Type of graph to plot. The default value is '*rates*'. Do not change!

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. plot_distest_perldoc:

.. _plot_distest_perldoc:

***************
plot_distest.pl
***************

Plot a figure with results from log-normal and negative binomial tests. 
Need to run :ref:`distribution_test.pl <distribution_test_perldoc>` and 
:ref:`grid_launcher_nb_test.pl <grid_launcher_nbtest_perldoc>` first to create 
test results. These test results should be in files:

	* WT_lnorm_lev_clean_test.dat
	* Snf2_lnorm_lev_clean_test.dat
	* WT_lnorm_lev_test.dat
	* Snf2_lnorm_lev_test.dat
	* WT_nb_deseq_lev_test.dat
	* Snf2_nb_deseq_lev_test.dat
	* WT_nb_lev_test.dat
	* Snf2_nb_lev_test.dat

::

  plot_distest.pl -norm=lev -psfile=dist_tests.ps
    
:param: `-norm` - Normalization used in distribution test scripts (see DESCRIPTION). This is used to identify the test result files. The default value is: *lev*.

:param: `-psfile` - Output postscript file.

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.
	
.. plot_fcp_perldoc:

.. _plot_fcp_perldoc:

***********
plot_fcp.pl
***********

Plot p-value versus fold-change. The input file for this script is a result 
from a DE tool, which should be a tab-delimited file containing a column 
with p-values and a column with log2 fold-change. A column with gene names 
allows for gene selection and highlighting::

	plot_fcp.pl -file=de_edger.txt -psfile=de_edger.ps     

:param: `-file` - A file with DE test results. As a minimum, it should contain a column with p-values and a column with log2 fold change. By default, the first column is a gene name, the second column is log2 fold change and the third column is the p-value. These can be changed by column options described below.

:param: `-psfile` - Name of the output postscript file.

:param: `-colg` - Column with gene name (identifier). The default value is 1.

:param: `-colfc` - Column with log2 fold change. The default value is 2.

:param: `-colp` - Column with p-value. The default value is 3.

:param: `-fcsign=[-1|1]` - Sign of log2 fold change. If the tool reports fold change the wrong way around, use -1 to reverse it.

:param: `-mcor=[none|bh|hb]` - Multiple test correction to apply in the plot.

  * *none* - the default value under the assumption that the tool/test output is already corrected
  * *bh* - Benjamini-Hochberg
  * *hb* - Holm-Bonferroni

:param: `-all` - If used, fold-change/p-value plots will be created for all tools in the file specified by *-testinfofile*.

:param: `-testinfofile` - A file with DE tools metadata. This is a text file containing records for each test/tool. Each record has a following format

::

	TEST deseq
	name = DEseq
	repfile = /cluster/gjb_lab/cdr/GRNAseq/analysis/bootstraps/DEseq/deseq_%dreps_defaultnorm.tsv
	sp_repfile = /cluster/gjb_lab/cdr/GRNAseq/analysis/bootstraps_with_spike_ins/deseq/deseq_%dreps_package_default.tsv
	fullfile = wrap_deseq.txt
	adjustment = bh
	FCsig = -1

*name* is the name of the test to be displayed. *repflie* is a format for 
bootstrap results (an FCP file) for all replicates (%d stands for replicate 
number). *sp_repfile* is for spike-in bootstrap results. *fullfile* is the 
result for full clean replicate set. *adjustment* is the multiple test 
adjustment used to obtain p-values. *FCsig* is an optional sign correction 
for fold-change. FCP files are tab-delimited with (at least) three columns: 
gene id, log2 fold change and p-value.

:param: `-sellist` - An optional comma-delimited list of files, containing gene names to highlight (in the first column). Genes from each file will be highlighted in a different colour.

:param: `-genlist` - An optional file with gene names (first column) to be displayed. Other genes will be ignored.

:param: `-fclines` - Comma-delimited list of fold-change values (not logarithms!) to be shown as vertical lines.

:param: `-title` - Title of the plot.

:param: `-alpha` - Significance level. The default value is 0.05. A dashed red line at the significance level will be shown in the plot and the number of the significantly DE genes, n_s, reported.

:param: `-lmin` - Minimum value (as a logarithm) for the vertical axis. If not specified, the script will calculate it from the data.

:param: `-fmax` - Minimum/maximum value (as a base-2 logarithm) for the horizontal axis. The default value is 3.

:param: `-label` - Optional label inside the plot (left-top corner).

.. plot_gene_examples_perldoc:

.. _plot_gene_examples_perldoc:

*********************
plot_gene_examples.pl
*********************

Plots a figure with expression distribution for several genes. Best-fitting 
log-normal and negative binomial distributions are shown::

	plot_gene_examples.pl -clean -psfile=figure.ps

:param: `-clean` - Use clean data replicates, as specified in *defs.dat* file.

:param: `-psfile` - Name of the postsript file to redirect output to.

.. plot_gene_pileup_perldoc:

.. _plot_gene_pileup_perldoc:

*******************
plot_gene_pileup.pl
*******************

Plots read pileup for a given gene/locus. It shows the mean pileup across 
replicates, plus/minus one standard deviation. Additionally, a replicate can 
be highlighted on top of the mean.

Pileup files need to be created using :ref:`build_pileup.pl <build_pileup_perldoc>`.
You also need to run :ref:`norms.pl <norms_perldoc>` first in order to created
normalization factor files::

	plot_gene_pileup.pl -gene=yhr215w -clean
	plot_gene_pileup.pl -locus=VIII:551800-553500 -clean
      
:param: `-gene` - Gene to plot. If not defined, *-locus* must be supplied instead.

:param: `-locus` - Locus to plot. Format is *chromosome:pos1-pos2*.

:param: `-clean` - If specified, only clean replicates will be used, as defined in *defs.dat*.

:param: `-norm` - Normalization to be used. Needs *<char>normfac.txt* files created with script :ref:`norms.pl <norms_perldoc>`. Default value is *deseq*.

:param: `-showrep` - Which replicate to show on top of mean and standard deviation of other replicates. The format is [r1],[r2] for condition 1 and 2. For example 21,6 will highlight replicate 21 in the first condition and 6 in the second one. You can also specify one replicate: 21, or ,6.

:param: `-psfile` - Output postscript file.

.. plot_ms_perldoc:

.. _plot_ms_perldoc:

**********
plot_ms.pl
**********

Plot mean, standard deviation and a few over things. This script can make 
lots of plots that we used to investigate data::

	plot_ms.pl
	plot_ms.pl -clean -withloess -psfile=figure.ps

:param: `-wthloess` - If specified, loess fit will be added to the plot.

:param: `-norm` - Expression data normalization. Input expression files contain raw counts. These are normalized on the fly, using one of the following methods:

  * *none* - no normalization
  * *lev* - levelling (equal read count) normalization (needs data preparation)
  * *deseq* - DESeq normalization (default)
  * *tmm* - trimmed mean of M values normalization
  * *fpkm* - approximate FPKM normalization (not identical to cuffdiff!)
  * *totcount* - total count
  * *totcountm* - total count with normalization factors scaled to the mean of 1

:param: `-clean` - Use clean data replicates, as specified in *defs.dat* file.

:param: `-psfile` - Name of the postsript file to redirect output to.

.. plot_powerstats_db_perldoc:

.. _plot_powerstats_db_perldoc:

*********************
plot_powerstats_db.pl
*********************

Plot various results::

	plot_powerstats_db.pl -test=edger     
    
:param: `-graph=[over|sim|cn|ns|sp]` - Which graph to plot. Only 'over' and 'sp' are working.

	* over - overview plot, used in the paper (default)
	* sp - significance propotion for three replicate numbers

:param: `-testinfofile` - A file with DE tools metadata. See :ref:`plot_fcp.pl <plot_fcp_perldoc>` for details. Default name is *de_tests.txt*.

:param: `-powerdir` - Directory where the results from :ref:`make_powerstats_db.pl <make_powerstats_db_perldoc>` are stored. The default value is *powerstats_db_ref*.

:param: `-test` - Name of the test/tool to plot results for. As defined in the metadata file.

:param: `-minnerp` - Minumum number of replicates to include in the plots. Default is 2.

:param: `-maxfc` - Maximum value of log2-fold-change axis in the plots. The default value is 2.  

:param: `-cs` - Larger character size in plots. The default value is 0.45. 

:param: `-cs2` - Smaller character size in plots. The default value is 0.3.

:param: `-psfile` - Name of the output postscript file.

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. plot_spikein_groups_perldoc:

.. _plot_spikein_groups_perldoc:

**********************
plot_spikein_groups.pl
**********************

::

	plot_fcp.pl -file=de_test_results.csv     
    
:param: `-file` - A file with DE test results. As a minimum, it should contain a column with p-values and a column with log2 fold change. By default, the first column is a gene name, the second column is log2 fold change and the third column is the p-value. These can be changed by column options described below.

:param: `-colg` - Column with gene name (identifier). The default value is 1.

:param: `-colfc>=I<number` - Column with log2 fold change. The default value is 2.

:param: `-colp>=I<number` - Column with p-value. The default value is 3.

:param: `-fcsign=[-1|1]` - Sign of log2 fold change. If the tool reports fold change the worng way around, use -1 to reverse it.

:param: `-all` - If used, fold-chage/p-value plots will be created for all tools in the file specified by *-testinfofile*.

:param: `-testinfofile` - A file with DE tools metadata. This is a text file contating records for each test/tool. Each record has a following format.

::

	TEST deseq
	name = DEseq
	repfile = /cluster/gjb_lab/cdr/GRNAseq/analysis/bootstraps/DEseq/deseq_%dreps_defaultnorm.tsv
	sp_repfile = /cluster/gjb_lab/cdr/GRNAseq/analysis/bootstraps_with_spike_ins/deseq/deseq_%dreps_package_default.tsv
	fullfile = wrap_deseq.txt
	adjustment = bh
	FCsig = -1

*name* is the name of the test to be displayed. *repflie* is a format for 
bootstrap results (an FCP file) for all replicates (%d stands for replicate 
number). *sp_repfile* is for spike-in bootstrap results. *fullfile* is the 
result for full clean replicate set. *adjustment* is the multiple test 
adjustment used to obtain p-values. *FCsig* is an optional sign correction 
for fold-change. FCP files are tab-delimited with (at least) three columns: 
gene id, log2 fold change and p-value.

:param: `-sellist` - An optional comma-delimited list of files, containing gene names to highlight (in the first column). Genes from each file will be highlighed in a different colour.

:param: `-genlist` - An optional file with gene names (first column) to be displayed. Other genes will be ignored.

..	toctree::
	:maxdepth: 1
