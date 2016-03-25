.. DEtoolcomp:

.. _DEtoolcomp:

***************************************
Comparing Differential Expression Tools
***************************************

.. gene_ranking_sigprop_perldoc:

.. _gene_ranking_sigprop_perldoc:

***********************
gene_ranking_sigprop.pl
***********************

Creates a table of gene ranking according to DE tests. Genes from each DE test 
are ranked, a score for each gene is calculated and a table (sorted by the 
score) with all genes from all tools is produced. This script uses a different 
approach to :ref:`gene_ranking.pl <gene_ranking_perldoc>`.  It looks at a 
proportion of significant DE calls across bootstraps for each gene. This is 
calculated by :ref:`make_powerstats_db.pl <make_powerstats_db_perldoc>` and 
stored in files *(prefix)_sigprop.stats*. This scripts simply reads these files 
and collates information::

	gene_ranking.pl -testinfofile=de_tests.dat -outfile=gene_ranking.txt   

:param: `-dir` - Directory with results from :ref:`make_powerstats_db.pl <make_powerstats_db_perldoc>`. 

:param: `-testinfofile` - A file with DE tools metadata. See :ref:`plot_fcp.pl <plot_fcp_perldoc>` for details.

:param: `-outfile` - Name of the output file. The default name is *gene_ranking.csv*.

:param: `-genlist` - A file with a list of gene names (first column) to use. The default name is *genlist.tsv*.

:param: `-nrep` - Use the specific number of replicates (from bootstrap tests) to do ranking. If not specified, the full clean replicate set results (as specified in *defs.dat* file) will be used.

:param: `-allrep` - Loop through all replicates and produce one ranking file for each.

:param: `-pfc` - If specified, ranking is done on p-values and fold changes (if not specified, ranking is done only on p-values).

:param: `-sig` - If specified, significance will be used to rank genes as opposed to p-value/fold change.

:param: `-multicor=[none|BH|HB]` - Multiple test correction to apply to tests/tools that return uncorrected raw p-values. The default value is *BH*.

:param: `-exclude` - A comma delimited list of tools/tests to exclude from analysis.
		
.. gene_ranking_perldoc:

.. _gene_ranking_perldoc:

***************
gene_ranking.pl
***************

Creates a table of gene ranking according to DE tests. Genes from each DE test 
are ranked, a score for each gene is calculated and a table (sorted by the 
score) with all genes from all tools is produced. There are two ways of 
ranking genes from individual DE tools (controlled by *-sig* option). 
The first method ranks them by the increasing p-values. Where p-values are 
identical (e.g. equal zero), ranking is (optionally, see option *-pfc*) done 
by the decreasing absolute fold change. The score for a gene is then the rank 
product.

The second method looks only at significance of each gene, 
according to a given criterion. 0/1 is assigned for significant/non-significant
gene. The score is the sum of these 0/1s::

	gene_ranking.pl -testinfofile=de_tests.dat -outfile=gene_ranking.txt   

:param: `-tests` - Comma-delimited list of test names (as defined in the file specified by *-testinfofile*) to be included.
    
:param: `-testinfofile` - A file with DE tools metadata. See :ref:`plot_fcp.pl <plot_fcp_perldoc>` for details.

:param: `-outfile` - Name of the output file. The default name is *gene_ranking.csv*.

:param: `-genlist` - A file with a list of gene names (first column) to use. The default name is *genlist.tsv*.

:param: `-pfc` - If specified, ranking is done on p-values and fold changes (if not specified, ranking is done only on p-values).

:param: `-sig` - If specified, significance will be used to rank genes as opposed to p-value/fold change.

:param: `-multicor=[none|BH|HB]` - Multiple test correction to apply to tests/tools that return uncorrected raw p-values. The default value is C<BH>.

.. make_powerstats_db_perldoc:

.. _make_powerstats_db_perldoc:

*********************
make_powerstats_db.pl
*********************

Calculates power statistics (true positives, true negatives, false positives, 
false negatives) for each test/tool against a chosen "true standard". The 
power stats are calculated for a range of replicates, while the true standard 
is a test result for the full clean data set. You need to run all power 
bootstraps first. The results are stored in a directory to be used later to 
create comparative power plots.

This script uses data directly from db files, and counts true/false 
positives/negatives for each bootstrap separately, then combining these 
results and reporting mean and standard deviation. An alternative approach 
is used in :ref:`make_powerstats.pl <make_powerstats_perldoc>`, where median 
p-values across bootstraps are used.

This script is a wrapper around 
:ref:`one_bs_powerstats.pl <one_bs_powerstats_perldoc>`, and it runs 
it on the cluster::

	make_powerstats_db.pl -test=edger

:param: `-test` - Name of the test/tool to plot results for. As defined in the metadata file.

:param: `-testinfofile` - A file with DE tools metadata. See :ref:`plot_fcp.pl <plot_fcp_perldoc>` for details. Default name is *de_tests.txt*.

:param: `-script` - This should be the full path to the :ref:`one_bs_powerstats_db.pl <one_bs_powerstats_db_perldoc>` script that this wrapper runs.

:param: `-logdir` - Path to the directory where all intermediate files are stored. Default is *./powerstats_logs*.

:param: `-powerdir` - Path to the directory where results are saved. Default value is *powerstats_db_ref.*

:param: `-maxn` - Maximum number of replicates to be used. The default value is 40.

:param: `-queue` - Grid engine queue name (or names, comma-separated).

:param: `-genlist` - A file with the list of genes. Created by :ref:`make_gene_list.pl <make_gene_list_perldoc>`. The default value is *genlist.tsv*.

:param: `-reffcfile` - A file with reference log2-fold-changes (in the second column) and gene names (first column). These fold changes will be used for consistency. This file can be created, e.g., by the following command. :command:`simple_de_test.pl -test=mw -clean -outfile=WT_Snf2_deseq_clean_mw_test.tsv`.

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.

.. one_bs_powerstats_perldoc:

.. _one_bs_powerstats_perldoc:

********************
one_bs_powerstats.pl
********************

To be used by :ref:`make_powerstats_db.pl <make_powerstats_db_perldoc>`. Not 
to be run on its own.

=head1 OPTIONS

=over 5

:param: `-genlist` - File with list of genes.

:param: `-dbfile` - Input sqlite .db file.

:param: `-totfile -fcfile1, -fcfile2, -rocfile, -truefile, -nsigfile, -sigpropfile` - Various temporary files to communicate with the parent script.

:param: `-multicor` - Multiple test correction.

:param: `-truecor` - Multiple test correction used in the 'true' file.

:param: `-alpha` - Significance limit. Default value is 0.05.

:param: `-fcsig` - Sign of log-fold-change.

:param: `-reffcfil` - The same as in the parent script.

:param: `-colp` - Column with p-values.

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.

..	toctree::
	:maxdepth: 1


