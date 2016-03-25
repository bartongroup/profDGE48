.. DEtools:
.. _DEtools:

************************************
Differential Expression Tool Scripts
************************************

These are the `R <https://www.r-project.org/>`_ scripts that actually run each 
Differential Expression (DE) call
on each of the bootstrap runs instances by :ref:`generic_wrapper.py <generic_wrapper_autodoc>`.
These scripts all have a standard interface, taking in four standard positional 
arguments that :ref:`generic_wrapper's R interface <genwrap_Rinterface>` provides.

The arguments are:

:param: `exprsFile` - A tab-delimited expression counts file. 
:param: `groupsFile` - A file containging the pheno, or grouping, information for the data
:param: `outFile` - The output filename for the called DE data
:param: `normalise` - Normalise the data? If this is set to *package:default'* the default normalization for the tool will be used.

Most of the tools are available in `Bioconductor <https://www.bioconductor.org/>`_. The analysis presented in `Gierlinski et. al. (2015) 
<http://bioinformatics.oxfordjournals.org/content/early/2015/07/22/bioinformatics.btv425.abstract>`_ 
and `Schurch et. al. (2015) <http://arxiv.org/abs/1505.02017>`_ used `R 3.2.2 <https://cran.r-project.org/src/base/R-3/R-3.2.2.tar.gz>`_
and `bioconductor 3.2 <https://www.bioconductor.org/packages/3.2/BiocViews.html>`_.

.. bayseq_Rdoc:

.. _bayseq_Rdoc:

********
baySeq.R
********

Run `baySeq v2.4.1 <https://www.bioconductor.org/packages/3.2/bioc/html/baySeq.html>`_ 
on a set of data::

	Rscript bayseq.R geneExprFile groupsFile outFile normalisation

.. cuffdiff_pq_perldoc:

.. _cuffdiff_pq_perldoc:

**************
cuffdiff_pq.pl
**************

Replace p-values (column 3) with q-values (column 4) in cuffdiff output

:param: `-infile` - Input file.

:param: `-outfile` - Output file.

.. degseq_Rdoc:

.. _degseq_Rdoc:

********
DEGSeq.R
********

Run `DEGSeq v1.24.0 <https://www.bioconductor.org/packages/3.2/bioc/html/DEGseq.html>`_ 
on a set of data::

	Rscript degseq.R geneExprFile groupsFile outFile normalisation

.. deseq_Rdoc:

.. _deseq_Rdoc:

*******
DESeq.R
*******

Run `DESeq v1.22.1 <https://www.bioconductor.org/packages/3.2/bioc/html/DESeq.html>`_ 
on a set of data::

	Rscript deseq.R geneExprFile groupsFile outFile normalisation

.. deseq2_Rdoc:

.. _deseq2_Rdoc:

********
DESeq2.R
********

Run `DESeq2 v1.10.1 <https://www.bioconductor.org/packages/3.2/bioc/html/DESeq2.html>`_ 
on a set of data::

	Rscript deseq2.R geneExprFile groupsFile outFile normalisation

.. ebseq_Rdoc:

.. _ebseq_Rdoc:

*******
EBSeq.R
*******

Run `EBSeq v1.10.0 <https://www.bioconductor.org/packages/3.2/bioc/html/EBSeq.html>`_ 
on a set of data::

	Rscript ebseq.R geneExprFile groupsFile outFile normalisation

.. edgeR_Rdoc:

.. _edgeR_Rdoc:

*******
edgeR.R
*******

Run `edgeR v3.12.0 <https://www.bioconductor.org/packages/2.11/bioc/html/edgeR.html>`_ 
on a set of data. This tool has two modes - exact (*edgeR.R*) and general
linear model (*edgeRglm.R*)::

	Rscript edgeR.R geneExprFile groupsFile outFile normalisation
	Rscript edgeRglm.R geneExprFile groupsFile outFile normalisation

.. limma_Rdoc:

.. _limma_Rdoc:

*******
limma.R
*******

Run `limma v3.26.8 <https://www.bioconductor.org/packages/3.2/bioc/html/limma.html>`_ 
on a set of data::

	Rscript limma.R geneExprFile groupsFile outFile normalisation

.. noiseq_Rdoc:

.. _noiseq_Rdoc:

********
noiseq.R
********

Run `NOISeq v2.14.1 <https://www.bioconductor.org/packages/3.2/bioc/html/NOISeq.html>`_ 
on a set of data::

	Rscript noiseq.R geneExprFile groupsFile outFile normalisation

.. poissonseq_Rdoc:

.. _poissonseq_Rdoc:

************
poissonseq.R
************

Run `PoissonSeq v1.1.2 <https://cran.r-project.org/web/packages/PoissonSeq/index.html>`_ 
on a set of data::

	Rscript poissonseq.R geneExprFile groupsFile outFile normalisation

.. run_cuffdiff_perldoc:

.. _run_cuffdiff_perldoc:

===============
run_cuffdiff.pl
===============
		
This script is a wrapper around `cuffdiff v2.1.1 <http://cole-trapnell-lab.github.io/cufflinks/releases/v2.1.1/>`_. 
It takes BAM files from provided directories, runs cuffdiff and converts its 
output to the required format::

	run_cuffdiff.pl -d=/cluster/gjb_lab/cdr/GRNAseq/mapping/genome_biol_reps -p=cuffdiff_results -g=Scerevisiae68_ERCC92.fasta -a=Saccharomyces_cerevisiae.EF4.68.gtf -o=cuffdiff_expr.tsv

:param: `-datapath|d` - Path to the data directory of the experiment. Each sub-directory in this will be treated as a separate condition in the experiment (with the name of the condition matching the name of the directory), and each .bam file in each directory is a replicate for that condition.

:param: `-outpath|p` - Path to the output directory, where cufflinks will write numerous output files.

:param: `-otufile|o` - The name (inc. path) of the output file from the wrapper.

:param: `-annotation|a` - Path to the .gff feature annotation file for the data. This file should match the feature you are counting the RNA-Seq expression for.

:param: `-genome|g` - Path to the FASTA file with reference genome. This is used by cuffdiff to run its bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates. You should have write access to the directory with this file as cufflinks will attempt to write an index file there. It might be a good idea to copy genome file (or a link) into your local directory.

:param: `-conditions|c` - Comma-separated list of conditions. (Default: *<WT,Snf2>*).

:param: `-threads` - Number of threads to be used by cuffdiff.

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.
		
	
.. samseq_Rdoc:

.. _samseq_Rdoc:

********
samseq.R
********

Run `SAMSeq v 2.0 <https://cran.r-project.org/web/packages/samr/>`_ on a set of data::

	Rscript samseq.R geneExprFile groupsFile outFile normalisation

.. ttest_equal_Rdoc:

.. _ttest_equal_Rdoc:

**************
t-test_equal.R
**************

Run a T-test on a set of data. Uses the `t.test() function <https://stat.ethz.ch/R-manual/R-patched/library/stats/html/t.test.html>`_ 
internal to R with var.equal set to *TRUE*, specifying equal variances in 
the data and uses a pooled variance estimate::

	Rscript t-test_equal.R geneExprFile groupsFile outFile normalisation

.. ttest_unequal_Rdoc:

.. _ttest_unequal_Rdoc:

****************
t-test_unequal.R
****************

Run a T-test on a set of data. Uses the `t.test() function <https://stat.ethz.ch/R-manual/R-patched/library/stats/html/t.test.html>`_ function internal to R with 
var.equal set to *FALSE*::

	Rscript t-test_unequal.R geneExprFile groupsFile outFile normalisation

..	toctree::
	:maxdepth: 1
