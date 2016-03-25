.. annotations:

.. _annotations:

***********
Annotations
***********

The annotations used to interpret the RNA-seq data for the project. Includes 
the ensembl release 68 annotations for Saccharomyces
cerevisiae and details of the ERCC artificial spike-in Mix kit from 
`ThermoFisher <https://www.thermofisher.com/order/catalog/product/4456740>`_.

.. ERCC_Controls_Analysis_tsvdoc:

.. _ERCC_Controls_Analysis_tsvdoc:

==========================
ERCC_Controls_Analysis.txt
==========================

Details of the 92 ERCC spike-in analysis values, including the concentrations 
of each spike-in in each mix.

.. ERCC_Controls_Annotation_tsvdoc:

.. _ERCC_Controls_Annotation_tsvdoc:

============================
ERCC_Controls_Annotation.txt
============================

Details of the 92 ERCC spike-ins, including ID, Genbank accessions and their 
sequence.

.. ERCC92_tsvdoc:

.. _ERCC92_tsvdoc:

==========
ERCC92.gtf
==========

Annotation information for the 92 ERCC spike-ins. (`gtf format 
<http://www.ensembl.org/info/website/upload/gff.html>`_)

.. get_gene_info_perldoc:

.. _get_gene_info_perldoc:

================
get_gene_info.pl
================

Get gene descriptions from Ensembl and store them in a tab-separated output file with three columns: gene id, gene name and description. The default name of the output file is C<gene_descriptions.tsv> in the current directory.::

	get_gene_info -genlist=WT_raw.tsv

:param: `-genlist` - Gene list, where first column contains gene names. Other columns are ignored. (default: `WT_raw.tsv`)

.. RNA_sample_details_tsvdoc:

.. _RNA_sample_details_tsvdoc:

======================
RNA_sample_details.tsv
======================

Details of the RNA samples used for the experiment. Includes Sample IDs, 
Saccharomyces cerevisiae strain, concentrations and RNA QC measurements

.. Saccharomyces_cerevisiae_with_spike_ins_tsvdoc:

.. _Saccharomyces_cerevisiae_with_spike_ins_tsvdoc:

==================================================
Saccharomyces_cerevisiae_with_spike_ins.EF4.68.gtf
==================================================

`Ensembl <http://www.ensembl.org/index.html>`_ release 68 Saccharomyces 
cerevisiae annotations with the 92 ERCC artificial spike-ins included 
(`gtf format <http://www.ensembl.org/info/website/upload/gff.html>`_)

.. Saccharomyces_cerevisiae_tsvdoc:

.. _Saccharomyces_cerevisiae_tsvdoc:

===================================
Saccharomyces_cerevisiae.EF4.68.gtf
===================================

`Ensembl <http://www.ensembl.org/index.html>`_ release 68 Saccharomyces 
cerevisiae annotations (`gtf format <http://www.ensembl.org/info/website/upload/gff.html>`_)

..	toctree::
	:maxdepth: 1
