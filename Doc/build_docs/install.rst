.. install:

.. _install:

************
Installation
************

None of the scripts in the :ref:`package <package_inventory>` require any
specific installation script or compilation script to be run. Just uncompress
them and away you go! 

That said, there are a considerable number of dependencies
for the full suite of scripts and they have only been tested with the specific
versions of the various software and :ref:`DE tools <DEtools>` listed. In 
addition many of them are optimized (or even require) access to a `DRMAA 
<http://www.drmaa.org/>`_ compatible linux cluster (preferably
an `SGE <http://gridscheduler.sourceforge.net/htmlman/manuals.html>`_ 
cluster, even more preferably the `University of Dundee, Life Sciences Compute 
Cluster <http://www.lifesci.dundee.ac.uk/services/lsc/services/cluster>`_).

The relevant environment will also need to be set up to allow access to each of
the installed dependencies. In particular, the environmental variables
PERL5LIB, PYTHONPATH, and PATH will need to be set correctly.

**Warning** some of the scripts (particularly the bash shell scripts) may 
contain hard-coded paths that will need to be modified for your install

=================
Software Versions
=================

  * `R v3.2.2 <https://cran.r-project.org/src/base/R-3/R-3.2.2.tar.gz>`_
  * perl `v5.10.1 <http://dev.perl.org/perl5/news/2008/perl-5.10.1.html>`_ / `v5.8.9 <http://dev.perl.org/perl5/news/2008/perl-5.8.9.html>`_
  * `python v2.6.6 <https://www.python.org/download/releases/2.6.6/>`_
  * `samtools v1.1.18 <http://sourceforge.net/projects/samtools/files/samtools/0.1.18/>`_
  * `TopHat v2.0.5 <http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.5.Linux_x86_64.tar.gz>`_
  * `BowTie2 v2.0.0-beta7 <http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.0.0-beta7/>`_
  * `Cufflinks v2.0.2 <http://cole-trapnell-lab.github.io/cufflinks/install/>`_

=============================
External Package Dependencies
=============================

  * python
  
    * `numpy <http://www.numpy.org/>`_
    * `scipy <http://www.scipy.org/>`_
    * `drmaa <http://drmaa-python.github.io/>`_
    * `pysqlite2 <https://pypi.python.org/pypi/pysqlite>`_
  
  * perl
  
    * `pdl <http://pdl.perl.org/>`_
    * `Ensembl release 68 API <http://jul2012.archive.ensembl.org/info/docs/api/api_installation.html>`_
    * `bioperl v1.2.3 <http://bioperl.org/DIST/old_releases/bioperl-1.2.3.tar.gz>`_
    * `Algorithm-cluster <http://search.cpan.org/~mdehoon/Algorithm-Cluster-1.52/perl/Cluster.pm>`_
    * `DBI <http://search.cpan.org/~timb/DBI-1.634/DBI.pm>`_
    * `Schedule-DRMAAc <http://search.cpan.org/~tharsch/Schedule-DRMAAc-0.81/Schedule_DRMAAc.pod>`_
    * `Time-HiRes <http://search.cpan.org/dist/Time-HiRes/HiRes.pm>`_
    * `Math-CDF <http://search.cpan.org/~callahan/Math-CDF-0.1/CDF.pm>`_
    * `IO-Compress <http://search.cpan.org/~pmqs/IO-Compress-2.068/>`_
    * `Thread <http://search.cpan.org/~jdhedden/threads-2.02/lib/threads.pm>`_
    * `Threads <http://search.cpan.org/~nwclark/perl-5.8.8/ext/threads/threads.pm>`_
    * `Carp <http://search.cpan.org/~rjbs/Carp-1.36/lib/Carp.pm>`_
    * `Exporter <http://search.cpan.org/~toddr/Exporter-5.72/lib/Exporter.pm>`_
    * `Devel-Size <http://search.cpan.org/~nwclark/Devel-Size-0.80/lib/Devel/Size.pm>`_
    * `Term-complete <http://search.cpan.org/~flora/Term-Complete-1.402/lib/Term/Complete.pm>`_
    * `Text-Tabs+Wrap <http://search.cpan.org/~muir/Text-Tabs%2BWrap-2013.0523/lib.old/Text/Wrap.pm>`_
  
  * R
  
    * `matrixStats <https://cran.r-project.org/web/packages/matrixStats/index.html>`_
    * `PoissonSeq <https://cran.r-project.org/web/packages/PoissonSeq/index.html>`_
    * `samr <https://cran.r-project.org/web/packages/samr/index.html>`_
    * `Bioconductor v2.11 <http://www.bioconductor.org/news/bioc_2_11_release/>`_

..	toctree::
	:maxdepth: 1

	getting_started

=====================================================
Dundee-specific Environment Variable settings & Paths
=====================================================

-----------------------
Environmental Variables
-----------------------

For :ref:`generic_wrapper.py <generic_wrapper_autodoc>` the Dundee-specific 
environmental variables are::

	export PYTHONPATH=*pathtocodebase*:/sw/opt/python/2.6.6/lib/python2.6/site-packages:/sw/opt/python/2.6.6/lib64/python2.6/site-packages
	export PERL5LIB=*pathtocodebase*/Modules:/sw/opt/ensembl-api/bioperl-live:/sw/opt/ensembl-api/68/ensembl/modules:/sw/lib/perl5.10.1/lib/perl:/sw/lib/perl5.10.1/lib/perl/x86_64-linux-thread-multi

For the downstram statisics and plotting scripts the Dundee-specific 
environmental variables are::

	export PERL5LIB=*pathtocodebase*/Modules:/homes/mgierlinski/perl/lib/perl5:/sw/perl/5.8.9/lib/perl5:/sw/lib/perl5/site_perl/5.8.9:/sw/lib/perl:/sw/lib/perl5/site_perl/5.8.9/x86_64-linux

-----
Paths
-----

These are specific paths that need to be provided to some of the scripts 
(notably :ref:`generic_wrapper.py <generic_wrapper_autodoc>`) for our setup.

  * samtools: /sw/opt/samtools-0.1.18/samtools
  * Rscript: /sw/opt/R/3.2.2/bin/Rscript
