.. Dists:

.. _Dists:

*************************************
Testing gene expression distributions
*************************************

.. distribution_test_perldoc:

.. _distribution_test_perldoc:

********************
distribution_test.pl
********************

Performs goodness-of-fit distribution test for gene expression data across 
biological replicates. For normal, log-normal and Poisson distributions tests 
are performed by this script. For across-lane Poisson and negative binomial, 
there are separate scripts, and this scritp only reads data. A figure is 
created with test results, showing distribution of p-values and p-value 
versus mean count plot::

	distribution_test.pl -dist=lnorm -clean -psfile=lnrom_test.ps
    
:param: `-dist` - Distribution to use for the test.

	* norm - normal
	* lnorm - log-normal
	* pois - Poisson
	* tpois - technical Poisson (between lanes); data are read from *all_poiss.dat* file, created with :ref:`poisson_reps.pl <poisson_reps_perldoc>`.
	* nb - negative binomial using the method described in `Meintanis (2005) <http://rivista-statistica.unibo.it/article/download/92/88>`_; you need to run :ref:`grid_launcher_nbtest.pl <grid_launcher_nbtest_perldoc>` first. Other NB methods do not work very well.

:param: `-lpoisfile` - File with lane goodness-of-fit results, created by :ref:`poisson_reps.pl <poisson_reps_perldoc>`.

:param: `-clean` - If defined, only clean replicates will be used, as defined in *defs.dat* file.

:param: `-norm=[none|totcount|totcountm|deseq|tmm|spikein]` - Normalization of gene expression data.

	* none - raw counts used - **not recommended**
	* totcount - total count divided by 1e6 - **not recommended**
	* totcountm - total count divided by its mean
	* deseq - DESeq normalization
	* tmm - TMM normalization
	* spikein - normalization to spike-ins (requires files with spike-in normalizing factors, *<condition>_spikein_normfac.txt* in countsdir directory, as defined in *defs.dat*)

:param: `-ptype=[dispersion|logratio]` - Type of Poisson test used. (Default: 'dispersion')

:param: `-minp` - Minimum p-value in the plot.

:param: `-out=[chauvenet|sigma]` - The method to reject outliers before doing test for normality or log-normality. If not specified, outliers are not rejected.

:param: `-sigma` - Sigma parameter for *sigma* method of rejecting outliers.

:param: `-nbstat=[m|a]` - Statistic to use with negative binomial test. 'm' - `Meintanis (2005) <http://rivista-statistica.unibo.it/article/download/92/88>`_, 'a' - `Anderson-Darling (1952) <http://projecteuclid.org/euclid.aoms/1177729437>`_. The latter one doesn't work very well. Default value is 'm'.

:param: `-report` - If specified a file will be created with mean and p-values, which can used for further analysis. 

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.

.. grid_launcher_nbtest_perldoc:

.. _grid_launcher_nbtest_perldoc:

***********************
grid_launcher_nbtest.pl
***********************

Runs negative binomial test (script :ref:`one_nb_test.pl <one_nb_test_perldoc>`) on the cluster::

	grid_launcher_nbtest.pl -script=./one_nb_test.pl -logdir=./nb_test -outfile=nbtest.txt

:param: `-script` - Location (full path and name) of the script :ref:`one_nb_test.pl <one_nb_test_perldoc>`.

:param: `-logdir` - Path to the directory where all intermediate files are stored.

:param: `-outfile` - Aggregated output file containing DE results.

:param: `-batchsize` - Number of genes in one batch. (Default: 30)

:param: `-ncrit` - Critical number of positive results to stop the simulation. See function *NBBootstrapTest* in module :ref:`Distribution.pm <Distribution_perlmod>`. (Default: 30)

:param: `-maxn` - Maximum number of simulations. (Default: 1e7)

:param: `-stat` - Statistic to use with Meintanis test. 'm' is for Meintanis, 'a' for Anderson-Darling. The default value is 'm'.

:param: `-options` - Any other options passed to the script.

:param: `-queue` - Cluster queue name(s), comma delimited.

:param: `-mem` - Memory requirement for each job. The default value is '4000M'.

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.
		
.. one_nb_test_perldoc:

.. _one_nb_test_perldoc:

**************
one_nb_test.pl
**************

This script is a part of :ref:`grid_launcher_nbtest.pl <grid_launcher_nbtest_perldoc>`.
It is not supposed to be run directly.

:param: `-batch` - Batch number.

:param: `-logdir` - Path to the directory where all intermediate files are stored.

:param: `-outfile` - Output file containing DE results.

:param: `-batchsize` - Number of genes in one batch. Default is 30.

:param: `-ncrit` - Critical number of positive results to stop the simulation. See function *NBBootstrapTest* in module :ref:`Distribution.pm <Distribution_perlmod>`. Default is 30.

:param: `-maxn` - Maximum number of simulations. Default is 1e7.

:param: `-stat` - Statistic to use with Meintanis test. 'm' is for Meintanis, 'a' for Anderson-Darling. The default value is 'm'.

:param: `-help` - Brief help.

:param: `-man` - Full manpage of program.

.. poisson_reps_perldoc:

.. _poisson_reps_perldoc:

***************
poisson_reps.pl
***************

Performs goodness-of-fit test for Poisson distribution across lanes. Requires 
gene expression files in a default directory (as defined in *defs.dat*) and a 
total read count file created with :ref:`count_fastq_reads_lane.pl <count_fastq_reads_lanes_perldoc>`::

	poisson_reps.pl  
    
:param: `-readcountfile` - The file with total read counts per lane, created with :ref:`count_fastq_reads_lane.pl <count_fastq_reads_lanes_perldoc>`. The default name is *all_readcount_lane.dat*.

:param: `-outfile` - The file with test results. The default name is *all_poiss.dat*. There are four columns written to this file: mean, standard deviation, p-value (not adjusted) and (chi-square distributed) T-statistic.

:param: `-ptype=[dispersion|logratio]` - Type of the test.

..	toctree::
	:maxdepth: 1


