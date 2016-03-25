.. normalisation:

.. _normalisation:

********************
Normalising the data
********************

.. levelling_lanes_perldoc:

.. _levelling_lanes_perldoc:

******************
levelling_lanes.pl
******************

Equal count (levelling) normalization of fastq files for lanes. Uses the 
grid engine.

  levelling_lanes.pl -levdir *path* 
    
:param: `-levdir` - Directory for levelled data.

:param: `-logdir` - Directory for log files.

:param: `-script` - Path and name of the script :ref:`one_levelling.pl <one_levelling_perldoc>`.

:param: `-queue` - Cluster queue name(s), comma delimited.

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.
		
.. levelling_perldoc:

.. _levelling_perldoc:

************
levelling.pl
************

Equal count (levelling) normalization of fastq files. No options other than 
those predefined in *defs.dat* so just run it. Read count files have to be 
prepared in advance using :ref:`count_fastq_reads.pl <count_fastq_reads_perldoc>`.

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. norms_perldoc:

.. _norms_perldoc:

********
norms.pl
********

Compares various normalization methods. Also, creates files with normalizing 
factors that are used by :ref:`plot_gene_pileup.pl <plot_gene_pileup>` and 
other scripts. There are no options.

.. one_levelling_perldoc:

.. _one_levelling_perldoc:

****************
one_levelling.pl
****************

Script used by :ref:`levelling_lanes.pl <levelling_lanes_perldoc>`. Not to be 
run on its own.
  
:param: `-infile` - Input file.

:param: `-outfile` - Output file.

:param: `-count` - Number of read counts.

:param: `-min` - Target number of reads.

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

..	toctree::
	:maxdepth: 1
