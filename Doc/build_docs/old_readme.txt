Scripts for data analysis
=========================

Data for this part of analysis are gene counts. These should be stored in two 
files (one per condition) with names WT_raw.tsv and Snf2_raw.tsv in countsdir 
directory. These can be created by "run_biological_alignemts.pl" and 
"group_by_gene.pl" script. All directories are defined in defs.dat file.


Preparing data
==============
combine_replicates.pl - combine individual replicate counts into multicolumn files, used by all other tools.
count_fastq_reads_lanes.pl - count reads in each condition, replicate and lane; creates a file required by Poisson test
count_fastq_reads.pl - count reads in each biological replicate, creates a file required by total count normalization
build_pileup.pl - build read pileup files for each chromosome (in Mapping Tools)
*pileup_exome.pl
pileup_var.pl - create read pileup statistic (required by bad_replicates.pl)
norms.pl - compare normalizations and store normalizing factors (required by plot_gene_pileup.pl)
levelling.pl - equal count normalization
levelling_lanes.pl - equal count normalization for lanes


Bad replicates
==============
compare_replicates.pl - various replicate comparisons (P1 Fig. 1)
bad_replicates.pl - create a figure that helps finding bad replicates (P1 Fig. 2)
*cluster_replicates.pl - cluster replicates using gene counts or read pileup profiles


Distribution tests
==================
poisson_reps.pl - Poisson test for lanes
*plot_lane_pvalues.pl - plot a distribution of p-value from Poisson test for lanes (requires poisson_reps.pl first)
grid_launcher_nbtest.pl - launches negative binomial tests on the cluster (requires one_nb_test.pl)
one_nb_test.pl - used by grid_launcher_nb_test.pl
distribution_test.pl - performs various distribution tests, including Poisson (requires output from poisson_reps.pl), log-normal and negative binomial (requires grid launcher first)
plot_distest.pl - plot a figure with results from log-normal and negative binomial tests (P1 Fig. 6)


Running basic DE tools
======================
simple_de_test.pl - a tool to run t-test, shrinkage t-test, log-ratio t-test or MW test 
grid_launcher_DE.pl - launches permutation and bootstrap tests on the cluster
one_bootstrap_ET93.pl - bootstrap test to be used with the grid launcher
one_permutest.pl - permutation test to be used with the grid launcher
run_cuffdiff.pl (in mapping tools) - wrapper around cuffdiff, required by grid_launcher_powertest_cuffdiff.pl


Sensitivity and power tests
===========================
generic_wrapper.py - run individual tools, bootstrap for increasing number of replicates
grid_launcher_powertest.pl - bootstraps for t-tests and MW test (requires one_powertest_simple.pl)
grid_launcher_powertest_cuffdiff.pl - bootstraps for cuffdiff (requires run_cuffdiff.pl)
make_cuffdiff_db.pl - create fold change/p-value databases from cuffdiff bootstrap results
combine_bs_pvalues.pl - calculate median p-values from bootstrap *.db files (requires one_bs_pvalue.pl)
one_bs_value.pl - used by combine_bs_values.pl


Tool comparison
===============
*compare_de_tests.pl - compare p-values or gene ranks from two DE tests
*compare_all_de_tests.pl - all comparisons in one plot
gene_ranking.pl - rank genes using DE tools results (default output 'gene_ranking.csv')
*gene_ranking_sigprop.pl - rank genes according to the proportion of significant bootstraps (cannot be used on full clean data set)
*plot_gene_cumul_ranking.pl - visualise gene ranking and compare tools
cluster_tools.pl - cluster DE tools according to their results
make_powerstats_db.pl - calculate power statistics (TP, TN, FP, FN) for each test/tool from bootstrap results, directly from db files (requires one_bs_powerstats.pl)
one_bs_powerstats.pl - used by make_powerstats_db.pl
plot_powerstats_db.pl - plot results from "make_powerstats_db.pl" (P2 Fig. 1, Figs.S2-S9)
compare_powerstats_db.pl - compare power test results from various tools (P2 Fig. 2)
compare_null.pl - compare results from null tests  - FPR (P2 Fig. 4)


Other plots
===========
plot_ms.pl - variance-mean and other plots (P1 Fig. 4)
*plot_hid.pl - Ratio-intensity diagram with selection of (significant) genes
*plot_one_gene.pl - plot read count distributions for a given gene
plot_gene_examples.pl - figure with count distribution from a few genes (P1 Fig. 5)
plot_fcp.pl - plot fold-change versus p-value (P2 Fig.S1)
plot_gene_pileup.pl - read pileup across replicates for a selected gene/locus (requires build_pileup.pl and norms.pl first) (P1 Fig. 3)

