#!/usr/bin/python

'''
***************************************************************************
The Great RNA-Seq Experiment Differential Expression Tool Bootstrap Wrapper
***************************************************************************

This file is a generic wrapper for running multiple bootstraps of
`R <https://www.r-project.org/>`_ Differential Expression (DE) algorithms.
Its written for
`python 2.6.4 <https://www.python.org/download/releases/2.6.4/>`_,
`R 3.2.2 <https://cran.r-project.org/src/base/R-3/R-3.2.2.tar.gz>`_ and
`bioconductor 3.2 <http://www.bioconductor.org/news/bioc_2_11_release/>`_.

This script interfaces with a suite of
`Rscript <https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html>`_
scripts - :ref:`one for each tool <DEtools>` - through a set of common input
parameters. This allows it to run each tool in a well-defined, consistent,
fashion for multiple random selections of the biological replicates in the
experiment. The script exposes the results through a standard `sqlite
<https://www.sqlite.org/>`_ database structure.

The idea of having a general wrapper for all the algorithms is to allow
us to run each of the DE algorithms in a consistent fashion, directly from the
the command line, and keep a careful track of exactly what tools were used at
every stage to get from the read data to the DE results.

The wrapper takes in a directory structure containing either:

  * bam files that are the result of aligning raw fastq's to a genome.
  * pre-generated (possibly from a previous run of this script!) `htseq-count
    <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_ output files.

This structure also defines the experimental conditions. If the input is a set
of bam files, it also requires a gtf/gff annotation file. It also requires
an output filename and log filename. The script is SGE cluster aware through the
`python drmaa library <http://code.google.com/p/drmaa-python/>`_ but can also
be run stand-alone. It takes a significant amount of time to process lots of
files and bootstraps though, so a cluster is recommended.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:version: 3.0

========
Overview
========

This section provides a gross overview of the operation of the
*generic_wrapper.py* script.

    1.  Parse command line arguments.
    2.  Gather script and program version information and relevant environmental
        variables and write these to the log file.
    3.  Parse experiment structure and get replicate file information.
    4.  Calculate or read in the gene readcounts. If the counts need to be
        generated from the bam files, the script uses the *perl* script
        :ref:`group_by_gene.pl <group_by_gene_perldoc>` which uses `htseq-count
        <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_ to
        aggregate the readcounts for each gene. If pre-generated counts are
        being used they are just read directly from the files.
    5.  Aggregate all the replicate gene expression data into a single matrix
        and begin the bootstrapping process.
    6.  For each bootstrap, if required, calculate total count normalization
        factors and normalize the data. This total read count for the
        replicates selected in the bootstrap iteration is calculated by a
        cluster-aware routine that calls `samtools
        <http://samtools.sourceforge.net/>`_ with the *flagstat* option, and
        parses the output. **Note:** we don't use this for the, calculations
        in :ref:`Schurch et. al. (2015) <paper2>`. There we use the default
        normalization method for each tool to run the DE analysis. These
        normalizations are aplied individually as part of :ref:`each tools
        Rscript <DEtools>`.
    7.  Write the (normalized) data to a file that is straightforward for R to
        read, and write the experiment structure to a file for R too!
    8.  Run the specified R script that loads the expression data, performs
        internal tool-specific normalization (if required), then runs the
        DE tools on the data and finally the writes the results to a sensible
        format file. Help for each of R scripts, including some background for
        each DE algorithm (should) be accessible in the `Differential
        Expression Algorithms <DEtools>`_ section.
    9.  Read and annotate the DE expression results with gene information from
        `ensembl <http://www.ensembl.org/index.html>`_ using the *perl* script
        :ref:`add_gene_name_column.pl <add_gene_name_column_perldoc>` to get
        the relevant annotation information via the *ensembl perl* API.
    10. Store the results in a `sqlite <https://www.sqlite.org/>`_ database
        file with a standardized structure and cleanup the temporary directory
        specified with *--tmpdir*.

======================
Command-line Arguments
======================

**usage\:**
    generic_wrapper.pl
    :option:`-d|--datapath` *<path>*
    :option:`-a|--annotation` *<file>*
    :option:`-r|--scriptfile` *<file>*
    :option:`-o|--outfile` *<file>*
    :option:`-l|--logfile` *<file>*
    [:option:`-s|species` *<str>*]
    [:option: `--precounts`
    [:option:`-f|--feature` *<str>*]
    [:option:`-c|--count-method` *<str>*]
    [:option:`-n|--norm` *<str>*]
    [:option:`--gbgfile` *<path><filename>*]
    [:option:`--agncfile` *<path><filename>*]
    [:option:`--Rpath` *<path>*]
    [:option:`--tmpdir` *<path>*]
    [:option:`--keep-tmpdir`]
    [:option:`--nocluster`]
    [:option:`--clustq` *<str>*]
    [:option:`--samtoolspath` *<path>*]
    [:option:`--restart`]
    [:option:`--verbose`]
    [:option:`--version`]
    [:option:`--help`]

.. gw_required:

.. _gw_required:

-------------------
Required Parameters
-------------------

.. option:: -d <path>, --datapath <path>

    Path to the data directory of the experiment. Each sub-directory in this
    will be treated as a separate condition in the experiment (with the name
    of the condition matching the name of the directory). Each .bam file or
    pre-calculated gene count file (defined by :option:`--precounts`)
    in each directory is treated as a replicate for that condition. BAM file
    indexes for each file should have the same filename, with .bai post-pended.

.. option:: -a <file>, --annotation <file>

    Path to the `GFF <https://www.sanger.ac.uk/resources/software/gff/>`_
    feature annotation file for the data. This file should
    match the features you are counting the RNA-Seq expression for. This is
    required even if you are using pre-calculated gene count files.

.. option:: -r <file>, --scriptfile <file>

    Name and path of Rscript file to execute in order to perform differential
    gene expression analysis. See :ref:`Differential Expression Algorithms
    <DEtools>` section for more details and help with each tool.

.. option:: -o <filename>, --outfile <filename>

    The name (inc. path if different fromt he current direstory) of the output
    sqlite file from the wrapper.

.. option:: -l <filename>, --logfile <filename>

    The name (inc. path) of the log file from the wrapper.

------------------
Annotation options
------------------

.. option:: -s <str>, --species <str> (Default: "Saccharomyces cerevisiae")

    The species used to generate the RNA-Seq data. This is used by the
    `ensembl <http://www.ensembl.org/index.html>`_ API to annotate the features
    found, so it must be a recognizable name format for `ensembl
    <http://www.ensembl.org/index.html>`_ (`this list
    <http://www.ensembl.org/info/about/species.html>`_ is probably a good place
    to start). Strings with spaces in must be in quotes.

----------------------------
reads-to-gene-counts options
----------------------------

.. option:: --precounts

    If this option is given the script looks for pre-generated gene count
    files in :option:`--datapath` The files to look for are defined by the
    :option:`--precount_fileext` options. If found, the gene counts from these
    files are used and gene count summarization from .bam files is skipped. If
    no files of the type are found, the script will look for bam files and
    will do gene count summarization if they are found.

    If this option and :option:`--savecounts` are set, and no gene count files
    are found, gene count summarization from the bam files will be done and the
    resulting gene count files saved to the path specified by
    :option:`--datapath`.

.. option:: --precount_fileext <str> (Default: gbgout)

    The file extension of pre-generated gene count data to look for.

.. option:: --savecounts

    If this option is given and the code summarises gene count information from
    bam file alignments, the resulting gene count files from running the gene
    count summarization will be saved to the path specified with
    :option:`--datapath`.

.. option:: -f <str>, --feature <str> (Default: gene_id)

    Sets the level at which you count the RNA-Seq signal. Options are:

        1. gene_id
        2. transcript_id

.. option:: -c <str>, --count-method <str> (Default: union)

    Choose the method used by `htseq-count
    <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_ for counting
    the reads for each feature. Options are:

        1. union
        2. istrict
        3. inoempty

---------------------
Normalization options
---------------------

.. option:: -n <str>, --norm <str> (Default: none)

    Normalization method for the tool being wrapped. Options are:

        1. none
        2. package:*normtype* - pass normalization control to the Rscript.
           *normtype* can be any of the possible normalization methods available
           from the Rscript for the tool.
        3. pmrm - Force Per Million Reads Mapped normalization.
        4. deseq - Force deseq-style normalisation.


-----------------
Bootstrap options
-----------------

.. option:: -k <int>, --kreps <int> (Default: None)

    The number of replicates to randomly select for each bootstrap. Replicates
    listed with :option:`--excludelist` are excluded from selection. Replicates
    listed with :option:`--includelist` are always included in the selection.
    If this option is *None*, or is a number greater than remaining number of
    replicates after excluding those specified with :option:`--excludelist`,
    then all the replicates are used. If this option is a number less than the
    remaining number of replicates fter excluding those specified with
    :option:`--excludelist`, then this number of replicates are selected at
    random for each bootstrap iteration.

.. option:: -b <int>, --nbootstrap <int> (Default: None)

    The number of bootstrap iterations to run.

.. option:: -i <file>, --includelist <file> (Default: None)

    A filename (including path is different from the current directory) for a
    text file listing (one per line) the full names of the replicate files
    to be explicitely included in every bootstrap selection.

.. option:: -e <file>, --excludelist <file> (Default: None)

    A filename (including path is different from the current directory) for a
    text file listing (one per line) the full names of the replicate files
    to be explicitely excluded from every bootstrap selection. This can be
    used to specify bad replicates to exclude. An example file is given in
    :ref:`exclude_badreps.tsv <exclude_badreps_tsvdoc>`.

-------------------------------------
External Dependancies and Environment
-------------------------------------

.. option:: --gbgfile (Default: ./group_by_gene.pl)

    The full path to the :ref:`group_by_gene.pl <group_by_gene_perldoc>` script.
    Called run `htseq-count
    <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_.

.. option:: --agncfile (Default: ./add_gene_name_column.pl)

    The full path to the :ref:`add_gene_name_column.pl
    <add_gene_name_column_perldoc>` script. Called to annotate the output from
    :ref:`group_by_gene.pl <group_by_gene_perldoc>`.

.. option:: --Rpath ((Default: /sw/opt/R/2.15.1/bin/Rscript)

    Specify the path to the Rscript executable for the R installation to use.

.. option:: --perlpath (Default: /user/bin/perl)

    Specify the path to the *perl* executable for the perl installation to use.

.. option:: --samtoolspath (Default: /local/bin/samtools)

    Specify the path to the `samtools <http://samtools.sourceforge.net/>`_
    executable to use.

---------------
Cluster options
---------------

.. option:: -j <int>, --njobs <int> (Default: 400)

    The max number of parallel cluster jobs to spawn if :option:`--nocluster`
    is not set.

.. option:: --clustq <str> (Default: 64bit-pri.q)

    Define a specific cluster queue to use to run cluster jobs.

.. option:: --nocluster

    If this option is given, the cluster support within the script will be
    diabled and everything will be run sequentially on the current machine.

-------------
Other options
-------------

.. option:: --tmpdir <path> (Default: ./.DEtemp)

    The path of the temporary directory you with to use for storing
    intermediate files.

.. option:: --keep-tmpdir

    If the entire script executes sucessfully the path defined with
    :option:`--tmpdir` is usually cleaned up. With this option set the
    cleanup is disabled and the temp directory and files are retained.

.. option:: --restart

    Restart a previously failed run. With this option set, this script will
    attempt to load the previous run settings and any parts of the run that
    were completed successfully from the *restart.pkl* file storred in
    the path provided with :option:`--tmpdir`.

.. option:: --version

    Show program's version number and exit.

.. option:: -v, --verbose

    Turn on verbose logging.

.. option:: -h, --help

    Print a basic description of the tool and its options to STDOUT.

----------------------------------------------------------
Internally-Passed options - **DO NOT SET THESE MANUALLY**
----------------------------------------------------------

.. option:: --bootstrapmaster

    Defines this as a master script that initiates bootstrap runs. This option
    should never be set manually. It is set internally to control spawned
    bootstrap cluster runs of this script.

.. option:: --bootstrapslave

    Defines this as a slave bootstrap run belonging to a set initiated by a
    :option:`--bootstrapmaster`. This option should never be set manually. It
    is set internally to control spawned bootstrap cluster runs of this script.

.. genwrapout:

.. _genwrapout:

======
Output
======

The output from a sucessful run of this script is a single `sqlite
<https://www.sqlite.org/>`_ database file with the following structure:

    * **Table:** *features* - features the DE is calculated for.

        Col: *id*

        Col: *featureID*

        Col: *name*

        Col: *desc*

    * **Table:** *bslogs* - the output logs from each of the bootstrap runs.

        Col: *id*

        Col: *log*

    * **Table:** *bootstraps* - details of the bootstraps performed.

        Col: *id*

        Col: *bsinlogid*

        Col: *comments*

        Col: *deruntime_sec*

        Col: *logid* - Foreign Key to *bslogs.id*

    * **Table:** *DEresults* - this table contains a column for each column of
      data output by the tool. To the best of our ability these results are
      output from the tools Rscript with a defined, consistent, set of columns.
      For the data presented in :ref:`Schurch et. al. (2015) <paper2>` this
      results in the columns described here:

        Col: *id*

        Col: *featureid* - Foreign Key to *features.id*

        Col: *bsid* - Foreign Key to *bootstraps.id*

        Col: *log2foldchange*

        Col: *significance*

        Col: *Snf2*

        Col: *WT*

These database files are then used by other scripts here to analyse and plot
the bootstrap results.

====================
Example Command-line
====================

This is an example command line used for testing. It runs :ref:`limma
<limma_Rdoc>` on a set of data for the ensembl 68 Saccharomyces cerevisiae
annotation, keeping all temp files, in verbose mode without using a cluster
looking for pre-generateded counts (saving them if they are have to be
calculated by the script from .bam files)::

    ./generic_wrapper.py -d testdata -a ../Annotations/Saccharomyces_cerevisiae.EF4.68.gtf
    -r limma.R --keep-tmpdir -v --gbgfile group_by_gene.pl --agncfile add_gene_name_column.pl
    -o test01.db -l test01.log --tmpdir testTMP --samtoolspath /sw/opt/samtools-0.1.18/samtools
    --precounts --nocluster --savecounts

'''

# IF YOU ADD AN IMPORT HERE, MAKE SURE YOU ADD IT TO THE 'imported_modules'
#LIST OR IT WON'T GET INCLUDED IN VERBOSE LOGGING OUTPUT!
imported_modules = ["os", "sys", "re", "time", "warnings", "cPickle",
                    "gzip", "numpy", "scipy", "textwrap", "optparse",
                    "subprocess", "drmaa", "math", "sqlite3", "shutil"]

try:
    import os, os.path, sys, re, time, warnings, cPickle, gzip, numpy, scipy
    import textwrap, optparse, subprocess, drmaa, math, shutil, tempfile
    from numpy.lib.recfunctions import append_fields, rename_fields
    from optparse import OptionParser, OptionGroup
    from pysqlite2 import dbapi2 as sqlite3
    from scipy.stats import mstats as mst
    from Modules.housekeeping import custom_formatwarning
    from Modules.housekeeping import parse_required_options
    from Modules.housekeeping import parse_allowed_values
    from Modules.housekeeping import parse_option_type
    from Modules.housekeeping import timeStr, createNewTempdir
    from Modules.housekeeping import write_log_header
except ImportError as e:
    print "\nGeneric Wrapper.py: Something went wrong when importing the " \
          "required modules. It's probably something to do with you not " \
          "having the right bits in your PYTHONPATH environmental variable.\n"
    print "Please make sure the path to the root of the checkout of in your " \
          "PYTHONPATH environmental variable.\n"
    print "If you're in Dundee, make sure you have path to the Modules " \
          "directory and /sw/lib/python2.6/site-packages in your PYTHONPATH. " \
          "Otherwise, please make sure your PYTHONPATH includes locations" \
          "for the following modules:\n\t%s\n" % ", ".join(imported_modules)
    print "Your PYTHONPATH is:"
    for line in sys.path:
        if line!="":
            print "\t%s" % line.strip()

    raise

version_string="3.0"

# force sequential output to the screen, rather than buffered output.
def printf(fmt, *args): sys.stdout.write(fmt % args)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

# override warning messages with a custom format
warnings.formatwarning = custom_formatwarning

def rep_subselect(experiment_structure, include, options, data=None):

    ''' subselects replicates from the full filename list

    :param dict experiment_structure: dictionary describing the directory, and
                                      hence the experiment, structure.

    :param list include: list of files to include in the subselection

    :param options: a valid :py:class:`optparse.Values` instance.

    :param data: a :py:class:`numpy.ndarray` structured array instance output
                 from :func:`read_expression_data`.

    :return: a tuple containing:

        * a new experiemnt structure with just the selected replicates in
        * a new replicate file list
        * either *None* or, if data has been provided, a new subset of data for
          just the selected replicates.

    The *options* instance must contain the keywords:

        * *kreps* - number of replicates to select for each botostrap iteration.
        * *log* - filename of the the log file to write to.

    Forces the inclusion of replicate filenames given in *include*. Called by
    :func:`parse_experiment_data`.

    '''

    # open log
    log = open(options.log,"a",0)
    log.write("\n %s: subselecting replicates..." % timeStr())

    sub_experiment_structure={}
    sub_replicate_file_list=[]
    sub_replicate_filename_list=[]
    for condition in experiment_structure:
        if len(experiment_structure[condition])<int(options.kreps):
            log.write("\n %s: Warning: Fewer replicates available in " \
                      "condition %s than requested in the --kreps option " \
                      "(%i). Returning all possible replicates." \
                      "" % (timeStr(),
                            condition,
                            int(options.kreps)
                            )
                      )
            sub_experiment_structure = experiment_structure
            sub_replicate_file_list = experiment_structure[condition]
        else:
            # first include all possible replicates that are in the include
            # list for this replicate
            these_increps=[]
            other_reps=[]
            for file in experiment_structure[condition]:
                if file in include or os.path.basename(file) in include:
                    these_increps.append(file)
                else:
                    other_reps.append(file)

            # check to see if we need to add more reps to this list
            if len(these_increps)>int(options.kreps):
                log.write("\n %s: Warning: condition %s contains more " \
                          "replicates specified in the include list " \
                          "than the number requested in --kreps (%i). " \
                          "Using all include replicates." \
                          "" % (timeStr(),
                                condition,
                                int(options.kreps)
                                )
                          )
            elif len(these_increps)<int(options.kreps):
                # randomly select additional repplicates
                other_reps = numpy.array(other_reps)
                numpy.random.shuffle(other_reps)
                i=0
                imax = int(options.kreps)-len(these_increps)
                while i<imax:
                    thisfile = os.path.basename(other_reps[i])
                    if thisfile not in sub_replicate_filename_list:
                        these_increps.append(other_reps[i])
                        sub_replicate_file_list.append(other_reps[i])
                        sub_replicate_filename_list.append(thisfile)
                    else:
                        imax+=1
                        if imax>len(other_reps):
                            msg = "\n %s: Warning: conditions contain " \
                                  "the same replicates (or, at least, " \
                                  "replicates with identical filenames) and " \
                                  "there are not enough unique replicates to " \
                                  "fill both conditions requested with " \
                                  "--kreps (%i)" \
                                  "" % (timeStr(), int(options.kreps))
                            log.write(msg)
                            raise RuntimeError(msg)

                    i+=1

            sub_experiment_structure[condition]=these_increps

    for key in sub_experiment_structure.keys():
        log.write("\n\t%s: selected %i replicates..." \
                  "" % (key, len(sub_experiment_structure[key])))
        for file in sorted(sub_experiment_structure[key]):
            log.write("\n\t\t%s" % (os.path.basename(file)))

    # subselect the expression data if provided
    sub_data = None
    if data is not None:
        log.write("\n %s: subselecting data...\n" % (timeStr()))
        replicate_index=[]
        for data_replicate in data[3]:
            data_replicate_bam = re.sub("\.gbgout_rawReads","",data_replicate)
            matched=False
            for file_replicate in sub_replicate_file_list:
                match = re.match(".+%s" % data_replicate_bam, file_replicate)
                if match:
                    matched=True
                    break
            replicate_index.append(matched)

        sub_exprs = data[0][:,numpy.where(replicate_index)[0]]
        sub_colids = list(numpy.array(data[3])[numpy.where(replicate_index)[0]])
        sub_data = (sub_exprs, data[1], data[2], sub_colids)

    log.close()
    return(sub_experiment_structure, sub_replicate_file_list, sub_data)

def runStructure(fileext, included, excluded, options):

    """ Parse the path for files & directory structure

    :param str fileext: extension of the files to look for

    :param list excluded: list of filenames to exclude

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * the experiment structure
        * a replicate file list
        * the list of forced included replicates

    The *options* instance must contain the *datapath* keyword giving the
    path to parse for data files to parse. The script walks the path
    characterising the directory structure and identifying files with a
    matching extension. It automatically excludes any files it finds that
    are on the *excluded* list. Called by :func:`parse_experiment_data`.

    """

    datapath_contents = os.listdir(options.datapath)

    experiment_structure = {}
    replicate_file_list=[]
    included_reps=[]
    excluded_reps=[]
    for thing in datapath_contents:
        if os.path.isdir("%s/%s" % (options.datapath,thing)):
            condition_files = os.listdir("%s/%s" % (options.datapath,thing))
            replicates = []
            for file in condition_files:
                if re.match(".+\.%s$" % fileext, file):
                    file_fullpath = "%s/%s/%s" % (options.datapath,thing,file)
                    poss_filenames = [file, file_fullpath,
                                      re.sub(".gbgout","",file),
                                      re.sub(".gbgout","",file_fullpath)]

                    if len(set(poss_filenames).intersection(set(excluded)))==0:
                        replicates.append("%s/%s/%s" % (options.datapath,
                                                        thing,
                                                        file
                                                        )
                                          )
                        replicate_file_list.append("%s/%s/%s" \
                                                   "" % (options.datapath,
                                                         thing,
                                                         file
                                                         )
                                                   )
                    else:
                        excluded_reps.append(file_fullpath)

            # mark files if they are in the includelist
            if len(included)>0:
                for file in replicates:
                    base_file = os.path.basename(file)
                    if base_file in included:
                        included_reps.append(base_file)

            experiment_structure[thing]=replicates

    return(experiment_structure, replicate_file_list, included_reps, excluded_reps)

def parse_experiment_data(options):

    ''' Gets condition and replicate information from the directory structure

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * a dictionary describing the experiment structure.
        * a list of all the replicates int eh experiemnt.
        * the list of forced included replicates.
        * a flag indicating whether or not to run gene count summarization
          on the files.

    The *options* instance must contain the keywords:

        * *datapath* - path of the data files to parse
        * *precounts* - toggle whether to look for pre-generated gene count
          files or use bam files.
        * *precount_fileext* - extension of pre-generated gene count files to
          look for.
        * *kreps* - number of replicates to select for each botostrap iteration.
        * *nbootstrap* - the number of bootstrap iterations to run.
        * *inc_reps* - text file containing a list (one per line) of the
          replicates that must be included in every bootstrap iteration.
        * *exc_reps* - text file containing a list (one per line) of the
          replicates that must be excluded from all bootstrap iterations.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.

    This subroutine forces the inclusion of files in the *options.inc_reps*
    file and excludes replicates on the *options.exc_reps*. If
    *options.nbootstrap* is none, but *options.kreps* is an integer,
    it will select *options.kreps* replicates once from the full list returned
    as *all_replicates*.

    '''

    # open log
    log = open(options.log,"a",0)

    include=[]
    exclude=[]

    if options.inc_reps is not None:
        log.write("\n %s: Parsing include list %s... " % (timeStr(),
                                                          options.inc_reps))
        file = open(options.inc_reps,"r")
        filedata = file.readlines()
        file.close()
        for file in filedata:
            include.append(file.strip())

    if options.exc_reps is not None:
        log.write("\n %s: Parsing exclude list %s... " % (timeStr(),
                                                          options.exc_reps))
        file = open(options.exc_reps,"r")
        filedata = file.readlines()
        file.close()
        for file in filedata:
            exclude.append(file.strip())

    do_genecounts=True
    datafile_ext = "bam"
    if options.precounts:
        do_genecounts=False
        datafile_ext = options.precount_fileext

    log.write("\n %s: Parsing %s structure looking for .%s files... " \
              "" % (timeStr(), options.datapath, datafile_ext))

    results = runStructure(datafile_ext, include, exclude, options)
    experiment_structure = results[0]
    replicate_file_list = results[1]
    included_reps = results[2]
    excluded_reps = results[3]

    # default to bam files if no matching pre-generated files are found
    if len(replicate_file_list)==0:
        if datafile_ext=="bam":
            raise IOError("No appropriate datafiles found.")
        else:
            log.write("\n %s: No .%s files found, trying .bam files... " \
                      "" % (timeStr(), datafile_ext))
            do_genecounts=True

            results = runStructure("bam", exclude, options)
            experiment_structure = results[0]
            replicate_file_list = results[1]
            included_reps = results[2]
            excluded_reps = results[3]

            if len(replicate_file_list)==0:
                raise IOError("No appropriate datafiles found.")

    log.write("found %i conditions:\n" % len(experiment_structure.keys()))

    all_files=[]
    for key in experiment_structure.keys():
        log.write("\t%s: found %i replicates...\n" \
                  "" % (key, len(experiment_structure[key])))
        for file in experiment_structure[key]:
            log.write("\t\t%s\n" % (os.path.basename(file)))
            if os.path.basename(file) not in all_files:
                all_files.append(os.path.basename(file))
            else:
                msg = "\n %s: Warning: the specified conditions contain " \
                      "replicates with identical filenames. Unless you are " \
                      "explicitely testing tool false positive rates (FPRs), " \
                      "this is not a good thing!! Unless you are testing " \
                      "FPRs, please rename the files in your conditions and " \
                      "ensure that the filenames are totally unique." \
                      "" % timeStr()
                log.write(msg)
                warnings.warn(msg)

    # check all include files are found
    if len(include)>0:
        if len(included_reps)!=len(include):
            log.write("\n %s: Error: Not all the replicates specified in the " \
                      "include list can be found in the directory structure." \
                      "Either one or more of the specified files don't exist," \
                      " or one or more files occur in both the exclude and " \
                      "include lists.\n\tFile\tFound?\n\t====\t======\n" \
                      "" % timeStr())
            for file in include:
                found = "no"
                if file in included_reps:
                    found="yes"
                log.write("\t%s\t%s\n" % (file, found))
            sys.exit()
    
    
    if len(exclude)>0:
        log.write("Excluded reps:\n")
        for exrep in excluded_reps:
            log.write("\t%s\n" % exrep)
    
    # subselect kreps from each condition, if nbootstrap is None
    if options.kreps is not None and options.nbootstrap is None:
        selection = rep_subselect(experiment_structure,
                                  include, options)
        experiment_structure = selection[0]
        replicate_file_list = selection[1]

    log.close()

    return(experiment_structure, replicate_file_list, include, do_genecounts)

def run_gbg_cluster(c_session, bam_file, options):

    ''' Instance a cluster run of gene count summarization for a bam file

    :param c_session: a :py:class:drmaa.Session: instance for interacting with
                      the cluster.

    :param str bam_file: the filename of the bam to summarise

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * the job ID number of the cluster job
        * the name of the file containing the resulting gene count summarization

    The *options* instance must contain the keywords:

        * *feature_annot* - the path to the `gff
          <https://www.sanger.ac.uk/resources/software/gff/>`_ file to use for
          gene summarization.
        * *cmeth* - the method for `htseq-count
          <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_ to use
          to summarize the read coutn data.
        * *feature* - the feature type in the annotation to summarise the read
          count data for.
        * *gbgpath* - the path to the :ref:`group_by_gene.pl
          <group_by_gene_perldoc>` script.
        * *perlpath* - the path to the perl executable to use to run
          :ref:`group_by_gene.pl <group_by_gene_perldoc>`.
        * *clustq* - the name of the cluster queue to submit jobs to.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.

    Instances a cluster run of :ref:`group_by_gene.pl <group_by_gene_perldoc>`
    for the given bamfile. Note that this function includes SGE-specific
    job submission parameters (in fact, they're quite possibly specific to our
    individual cluster!). Briefly, they are:

        * *-cwd* - specify running the job in the current working directory
        * *-l local_free=20G* - require at least 20GB of free disk local to the
          cluster node
        * *-l ram=3G* - require at least 3GB of RAM on the node.
        * *-clear* - clear any existing cluster queue settings.
        * *-q 'options.clustq'* - specify the queue to submit the job to.

    '''

    # open log
    log = open(options.log,"a",0)
    log.write("\n\tProcessing file %s..." % bam_file)

    outfile = "%s/%s.gbgout" % (options.tmpdir,
                                os.path.basename(bam_file))

    c_job = c_session.createJobTemplate()

    c_job.remoteCommand = options.perlpath
    c_job.args = [options.gbgpath, "--bam", bam_file, "--gtf",
                  options.feature_annot, "--method", options.cmeth,
                  "--feature", options.feature, "--out", outfile]

    if options.verbose:
        c_job.args.append("--debug")
        log.write("\n\t\tdrmaa command line:\n\n%s %s\n" \
              "" % (options.perlpath," ".join(c_job.args)))

    c_job.outputPath = ":%s" % options.tmpdir
    c_job.errorPath = ":%s" % options.tmpdir

    # pass current working directory (not that this is needed really, but hey!)
    c_job.nativeSpecification = "-cwd"

    # check for local disk space on the target node
    c_job.nativeSpecification = "%s -l local_free=20G -l ram=3G" \
                                "" % c_job.nativeSpecification

    # add support for different cluster queue specifications
    c_job.nativeSpecification = "-clear -q '%s' %s" \
                                "" % (options.clustq, c_job.nativeSpecification)

    if options.verbose:
        log.write("\n\t\tdrmaa output intermediates written to: %s" \
                  "" % options.tmpdir)

    c_job.jobEnvironment = os.environ
    jobid = c_session.runJob(c_job)

    log.write("\n\t\tJob submitted with id: %s\n" % jobid)

    log.close()

    return(jobid, outfile)

def run_gbg(bam_file, options):

    ''' Instance a non-cluster run of gene count summarization for a bam file

    :param str bam_file: the filename of the bam to summarise

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * *None*
        * the name of the file containing the resulting gene count summarization

    The *options* instance must contain the keywords:

        * *feature_annot* - the path to the `gff
          <https://www.sanger.ac.uk/resources/software/gff/>`_ file to use for
          gene summarization.
        * *cmeth* - the method for `htseq-count
          <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_ to use
          to summarize the read coutn data.
        * *feature* - the feature type in the annotation to summarise the read
          count data for.
        * *gbgpath* - the path to the :ref:`group_by_gene.pl
          <group_by_gene_perldoc>` script.
        * *perlpath* - the path to the perl executable to use to run
          :ref:`group_by_gene.pl <group_by_gene_perldoc>`.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.

    The script uses :py:class:`subprocess.Popen` to call an external
    run of :ref:`group_by_gene.pl <group_by_gene_perldoc>`.

    '''

    # open log
    log = open(options.log,"a",0)
    log.write("\tProcessing file %s...\n" % bam_file)

    args = [options.perlpath, options.gbgpath, "--bam", bam_file, "--gtf",
            options.feature_annot, "--method", options.cmeth, "--feature",
            options.feature, "--out",
            "%s/%s.gbgout" % (options.tmpdir,
                              os.path.basename(bam_file))]

    if options.verbose:
        log.write("\t\tcommand line:\n\t\t%s\n" \
              "" % " ".join(args))
        log.write("\t\tJob executing...\n")

    gbg_call = subprocess.Popen(args, stdout=subprocess.PIPE)
    gbg_out = gbg_call.communicate()[0]
    log.write("\t\tcompleted:\n%s\n" % gbg_out)

    log.close()

    return(None, "%s/%s.gbgout" % (options.tmpdir, os.path.basename(bam_file)))

def run_agnc(agnc_file, options, coords=True, desc=True):

    ''' Instance a non-cluster run of gene count summarization annotation

    :param str agnc_file: the filename of the gene count summarization file to
                          annotate

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * *None*
        * the name of the file containing the resulting annotated gene count
          summarization

    The *options* instance must contain the keywords:

        * *species* - the ensembl species name.
        * *agncpath* - the path to the :ref:`add_gene_name_column.pl
          <add_gene_name_column_perldoc>` script.
        * *perlpath* - the path to the perl executable to use to run
          :ref:`add_gene_name_column.pl <add_gene_name_column_perldoc>`.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.

    The script uses :py:class:`subprocess.Popen` to call an external
    run of :ref:`add_gene_name_column.pl <add_gene_name_column_perldoc>`.
    :ref:`add_gene_name_column.pl <add_gene_name_column_perldoc>` uses the
    `perl ensembl API <http://www.ensembl.org/info/docs/api/index.html>`_
    to pull species annotation information directly from ensembl. *species*
    must be a name recognised by ensembl (`this list
    <http://www.ensembl.org/info/about/species.html>`_ is probably a good place
    to start). Strings with spaces in must be in quotes.

    '''


    # open log
    log = open(options.log,"a",0)
    log.write("\n\tProcessing file %s..." % agnc_file)

    args = [options.perlpath, options.agncpath, "--in", agnc_file,
            "--species", options.species,
            "--out", "%s/%s.agncout" % (options.tmpdir,
                                     os.path.basename(agnc_file))
            ]

    if coords:
        args.append("--coords")

    if desc:
        args.append("--desc")

    if options.verbose:
        log.write("\n\t\tcommand line:\n\t\t%s" \
              "" % " ".join(args))
        log.write("\n\t\tJob executing...")

    agnc_call = subprocess.Popen(args, stdout=subprocess.PIPE)
    agnc_out = agnc_call.communicate()[0]
    log.write("completed.")

    log.close()

    return(None, "%s/%s.agncout" % (options.tmpdir,
                                    os.path.basename(agnc_file)))

def monitorDrmaaJobs(c_session, joblist, filelist, prefix, options, delo=True,
                     dele=True):

    """ monitor a list of cluster jobs that have been submitted via drmaa

    :param c_session: a :py:class:drmaa.Session: instance for interacting with
                      the cluster.

    :param list joblist: a list of cluster job IDs to monitor.

    :param list filelist: a list of the files to expect from sucessful
                          completion of each monitored job.

    :param str prefix: the cluster output and error file prefix to look for.

    :param options: a valid :py:class:`optparse.Values` instance.

    :param bool delo: toggle deletion of the cluster output file on job
                      completion.

    :param bool dele: toggle deletion of the cluster error file on job
                      completion.

    :return: a tuple containing:

        * A dictionary keyed on the job ID giving the completion status of each
          monitored job.
        * the number of jobs that completed successfully
        * the number of jobs that failed to complete successfully.
        * the input (possibly modified) :py:class:`optparse.Values` instance.

    The *options* instance must contain the keywords:

        * *log* - filename of the the log file to write to.
        * *keep_tmpdir* - flag to keep or delete the tmpdir at the end of the
          script run.
        * *verbose* - boolean toggle for verbose logging.

    At the moment *filelist* si limited to one file per job and so this list
    must have the same length and be in the same order as *joblist*. As well as
    monitoring the jobs as they run, and looking for the expected output files,
    the script also looks for the cluster output files. They are expected to
    have the forms:

        * *<prefix>.o<jobid>* - cluster output log.
        * *<prefix>.e<jobid>* - cluster error log.

    If any of the runs fail, *options.keep_tmpdir* will be overwritten to
    **True** to ensure all failing cluster log files are kept.

    """

    # open log
    log = open(options.log,"a",0)
    log.write("\n\t\tmonitoring jobs...")

    jobs_remaining = {}
    jobs_to_files = {}
    i=0
    while i<len(joblist):
        jobs_remaining[joblist[i]] = None
        jobs_to_files[joblist[i]] = filelist[i]
        i+=1

    finished=False
    exit_status = {}

    success = 0
    fail = 0
    while not finished:
        for jid in jobs_remaining:
            if jid not in exit_status.keys():
                status = c_session.jobStatus(jid)
                if jobs_remaining[jid] is None or jobs_remaining[jid]!=status:

                    ofilename = "%s.o%s" % (prefix, str(jid))
                    efilename = "%s.e%s" % (prefix, str(jid))

                    # cluster job status checking...
                    # Job has finished normally
                    if status==drmaa.JobState.DONE:

                        exit_status[jid]=status
                        success+=1

                        try:
                            sizecheck = os.path.getsize(jobs_to_files[jid])
                            if sizecheck>0:
                                try:
                                    if delo and dele:
                                        if options.verbose:
                                            log.write("\n %s: %s - completed sucessfully." \
                                                      "\n\t removing file: %s" % (timeStr(),
                                                                                  str(jid),
                                                                                  ofilename)
                                                      )
                                        os.remove(ofilename)
                                        if options.verbose:
                                            log.write("\n\t removing file: %s" % efilename)
                                        os.remove(efilename)
                                    elif delo and not dele:
                                        if options.verbose:
                                            log.write("\n %s: %s - completed sucessfully." \
                                                      "\n\t removing file: %s" % (timeStr(),
                                                                                  str(jid),
                                                                                  ofilename)
                                                      )
                                        os.remove(ofilename)
                                    elif dele and not delo:
                                        if options.verbose:
                                            log.write("\n %s: %s - completed sucessfully." \
                                                      "\n\t removing file: %s" % (timeStr(),
                                                                                  str(jid),
                                                                                  efilename)
                                                      )
                                        os.remove(efilename)
                                except:
                                    log.write("\n\tremoving files for sucessful job %s " \
                                              "failed." % str(jid)
                                              )
                            else:
                                # keep all tmpfiles regardless of settings and report error
                                options.keep_tmpdir = True
                                exit_status[jid] = "emptyfile"
                                log.write("\n %s: %s: Output file %s is empty." \
                                          " It might just be taking a while to " \
                                          "write the file but, I'll keep all " \
                                          "the temp files incase something went " \
                                          "wrong." % (timeStr(), str(jid), jobs_to_files[jid])
                                          )
                        except:
                            # keep all tmpfiles regardless of settings and report error
                            options.keep_tmpdir = True
                            exit_status[jid] = "nofile"
                            log.write("\n %s: Output file %s does not exist! Keeping all " \
                                      "temp files incase something went wrong." \
                                      "" % (timeStr(),jobs_to_files[jid])
                                      )
                    # Job finished, but terminated abnormally
                    elif status==drmaa.JobState.FAILED:
                        # keep all tmpfiles regardless of settings and report error
                        options.keep_tmpdir = True
                        fail+=1
                        exit_status[jid]=status
                        log.write("\n %s: %s - job failed. Files for this " \
                                  "job are: \n\t%s\n\t%s" % (timeStr(), jid,
                                                    ofilename, efilename)
                                  )
                    # Job is queued and active
                    elif status==drmaa.JobState.QUEUED_ACTIVE:
                        if options.verbose:
                            log.write("\n %s: %s - job queued and active..." \
                                      "" % (timeStr(),jid)
                                      )
                    # job is running
                    elif status==drmaa.JobState.RUNNING:
                        if options.verbose:
                            log.write("\n %s: %s - job running..." \
                                      "" % (timeStr(),jid)
                                      )
                    # job is queued and in system hold
                    elif status==drmaa.JobState.SYSTEM_ON_HOLD: #
                        if options.verbose:
                            log.write("\n %s: %s - job on hold (by system)..." \
                                      "" % (timeStr(),jid)
                                      )
                    # job is queued and in user hold
                    elif status==drmaa.JobState.USER_ON_HOLD:
                        if options.verbose:
                            log.write("\n %s: %s - job on hold (by user - " \
                                      "why??)..." % (timeStr(),jid)
                                      )
                    # job is system suspended
                    elif status==drmaa.JobState.SYSTEM_SUSPENDED:
                        if options.verbose:
                            log.write("\n %s: %s - job suspended (by system)" \
                                      "..." % (timeStr(),jid)
                                      )
                    # job is user suspended
                    elif status==drmaa.JobState.USER_SUSPENDED:
                        if options.verbose:
                            log.write("\n %s: %s - job suspended (by user - " \
                                      "why??)..." \
                                      "" % (timeStr(),jid)
                                      )
                    # job is queued and in user and system hold
                    elif status==drmaa.JobState.USER_SYSTEM_ON_HOLD:
                        if options.verbose:
                            log.write("\n %s: %s - job on hold (by system & " \
                                      "user - wtf?? seriously??)..." % (timeStr(),jid)
                                      )
                    # process status cannot be determined,
                    elif status==drmaa.JobState.UNDETERMINED:
                        # keep all tmpfiles regardless of settings and report error
                        options.keep_tmpdir = True
                        fail+=1
                        exit_status[jid]=status
                        log.write("\n %s: %s - job state undetermined, " \
                                  "assumed success, but keeping files just " \
                                  "in case! Files for this job are: " \
                                  "\n\t%s\n\t%s" % (timeStr(),
                                                    jid,
                                                    ofilename,
                                                    efilename)
                                  )
                    else:
                        sys.exit("Something went wrong with the drmaa/cluster " \
                                 "monitoring. Get help!")

                    jobs_remaining[jid]=status

        for key in exit_status:
            try:
                del jobs_remaining[key]
            except:
                pass

        if len(jobs_remaining.keys())==0:
            finished=True

    log.close()

    return(exit_status, success, fail, options)

def get_gene_counts(all_replicates, options):

    ''' Get gene count summarization for all bam files in the experiemnt

    :param list all_replicates: a list of all the replicates in the experiment
                                (2nd element in the tuple output from
                                :func:`parse_experiment_data`).

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a list of the resulting file containing the gene count
             summarisation from each bam file

    The *options* instance must contain the keywords:

        * *feature_annot* - the path to the `gff
          <https://www.sanger.ac.uk/resources/software/gff/>`_ file to use for
          gene summarization.
        * *cmeth* - the method for `htseq-count
          <http://www-huber.embl.de/users/anders/HTSeq/doc/count.html>`_ to use
          to summarize the read coutn data.
        * *feature* - the feature type in the annotation to summarise the read
          count data for.
        * *gbgpath* - the path to the :ref:`group_by_gene.pl
          <group_by_gene_perldoc>` script.
        * *perlpath* - the path to the perl executable to use to run
          :ref:`group_by_gene.pl <group_by_gene_perldoc>`.
        * *savecounts* - toggle saving the gene count summarization files.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.
        * *nocluster* - toggle running the summarization on the cluster, or
          locally.


    Calls either :func:`run_gbg` or :func:`run_gbg_cluster` depending on
    whether *options.nocluster* is *True* or *False*. If *options.nocluster* is
    *False* then then *options* must also contain the keyword:

        * *clustq* - the name of the cluster queue to submit jobs to.

    '''

    # run group_by_gene.pl on each of the bam files. Use cluster is available
    gbg_files={}

    if options.nocluster:
        # open log
        log = open(options.log,"a",0)
        log.write("\n %s: Cluster not enabled: Running group_by_gene.pl and " \
                  "add_gene_name_column.pl on each bam file in the " \
                  "experiment...\n" \
                  "" % timeStr())
        log.close()

        for bam_file in all_replicates:
            job, gbg_outfile = run_gbg(bam_file, options)
            gbg_files[gbg_outfile]=True

    else:
        # open log
        log = open(options.log,"a",0)
        log.write("\n %s: Cluster enabled: Running group_by_gene.pl and " \
                  "add_gene_name_column.pl on a cluster node in the %s " \
                  "queue for each bam file in the experiment...\n" \
                  "" % (timeStr(),options.clustq))
        log.close()

        # use drmaa to spread multiple group_by_gene.pl jobs across the
        # cluster.
        # open a drmaa
        c_session = drmaa.Session()
        c_session.initialize()

        # OK, so samtools seems to have an issue with our cluster and
        # occasionally, randomly dies producing no error, just a core dump
        # we hypothesis this is due to OOM problems when the cluster is busy
        # so we're going to check for any failed jobs here and force a resubmit.
        resub_trys=1
        resub_count=0

        all_replicate_save = all_replicates

        while len(all_replicates)>0:

            if resub_count>0:
                log = open(options.log,"a",0)
                log.write("\n %s: resubmitting %i %s jobs that failed " \
                          "inexplicably..." % (timeStr(),
                                               len(all_replicates),
                                               os.path.basename(options.gbgpath)
                                               )
                          )
                log.close()

            # loop through all the replicates and save only those ones with
            # unique filenames. Here we are assuming that if two files in
            # different conditions have the same filename, they are the same!
            # So don't do this if they are not.

            baselist = []
            process_files=[]
            for bam_file in all_replicates:
                base_filename = os.path.basename(bam_file)
                if base_filename not in baselist:
                    baselist.append(base_filename)
                    process_files.append(bam_file)

            # loop throught he files calling a new cluster job for each
            joblist = []
            job_bamfile={}
            for bam_file in process_files:
                jobid, gbg_outfile = run_gbg_cluster(c_session, bam_file,
                                                     options)
                joblist.append(jobid)
                job_bamfile[jobid]=bam_file
                gbg_files[gbg_outfile]=True

            log = open(options.log,"a",0)
            log.write("\n %s: %s jobs launched..." \
                      "" % (timeStr(), os.path.basename(options.gbgpath)))
            log.close()

            prefix = "%s/%s" % (options.tmpdir,
                                os.path.basename(options.perlpath))
            filelist = gbg_files.keys()
            exit_status, sucess, fail, options = monitorDrmaaJobs(c_session,
                                                                  joblist,
                                                                  filelist,
                                                                  prefix,
                                                                  options)

            # clearup job metadat and information
            for jobid in joblist:
                c_session.synchronize(jobid,
                                      drmaa.Session.TIMEOUT_WAIT_FOREVER,
                                      True)

            # reset all_replicates list
            all_replicates=[]

            # only resubmit a fixed number of times...
            if resub_count<resub_trys:
                # check for any jobs that died without
                # producing a file or producing an empty file
                for entry in exit_status:
                    status = exit_status[entry]
                    if status=="emptyfile" or status=="nofile":
                        # recheck, just in case it took a while to write the
                        # file:
                        try:
                            sizecheck = os.path.getsize(job_bamfile[entry])
                            if sizecheck==0:
                                all_replicates.append(job_bamfile[entry])
                        except:
                            all_replicates.append(job_bamfile[entry])

            resub_count+=1

        c_session.exit()
        all_replicates = all_replicate_save

    # copy the gene counts to the source directory if specified.
    if options.savecounts:
        log = open(options.log,"a",0)
        log.write("\n %s: Saving gene read counds to source datafile " \
                  "directory..." % timeStr())
        log.close()
        for filename in gbg_files.keys():
            if options.verbose:
                log = open(options.log,"a",0)
                log.write("\n %s: Searching for match for file %s." \
                          "" % (timeStr(),os.path.basename(filename)))
                log.close()
            for origin_file in all_replicates:
                log = open(options.log,"a",0)
                log.write("\n %s: matching against %s..." \
                          "" % (timeStr(),os.path.basename(origin_file)))
                log.close()
                if re.match("^%s.gbgout$" % os.path.basename(origin_file),
                            os.path.basename(filename)):
                    try:
                        log = open(options.log,"a",0)
                        log.write("\n %s: Saving file %s to " \
                                  "locations %s." \
                                  "" % (timeStr(),
                                        filename,
                                        os.path.dirname(origin_file)))
                        log.close()
                        shutil.copy2(filename, os.path.dirname(origin_file))
                    except IOError, e:
                        log = open(options.log,"a",0)
                        log.write("\n %s: Problem copying file %s to " \
                                  "locations %s. Its probably a write " \
                                  "permissionas issue!" \
                                  "" % (timeStr(),
                                        filename,
                                        os.path.dirname(origin_file)))
                        log.close()
                        raise(e)

    return(gbg_files.keys())

def read_expression_data(replicate_files, options, format="unstruct",
                         delimiter="\t"):

    ''' reads the gene count summarization data & concatenates to a numpy array

    :param list replicate_files: a list of file containing replicate gene count
                                 data (output from :func:`get_gene_counts`).

    :param options: a valid :py:class:`optparse.Values` instance.

    :param str format: the format of the output data from this function.
                       Options are:

                       * unstruct - return a tuple of information
                       * struct - return a :py:class:`numpy.ndarray` structured
                         array

    :param str delimiter: The delimiter to use for separating fields in the
                          read count summarization files.

    :return: The functions returns either a :py:class:`numpy.ndarray` structured
             array with feature, coordinate and description columns that
             contains the gene count summarization data for the
             replicates, or a tuple containing:

               * a :py:class:`numpy.ndarray` with the read count summarization
                 data.
               * a :py:class:`numpy.ndarray` containing the feature names the
                 read count summarization applies to (in the same order as the
                 rows of the data array).
               * *None* (don't ask).
               * a list of the samples the read counts were summarized for (in
                 the same order as the columns of the data array).

             The default is *unstruct*.

    The *options* instance must contain the *log* keyword providing the
    filename of the the log file to write to.

    '''

    # open log
    log = open(options.log,"a",0)
    log.write("\n %s: Concatenating data in all files...\n" \
                  "" % timeStr())

    expressiondata = None
    file_structure = {'names':('featureID', 'rawReads'),
                      'formats':('S50', 'int')}

    datacolumns=[]
    for file in replicate_files:
        filename = os.path.basename(file)
        log.write("\tprocessing %s...\n" % filename)
        filedata = numpy.genfromtxt(file,
                                    dtype=file_structure,
                                    delimiter=delimiter,
                                    comments="$(*&$(*&$(*" # allow strings w. #
                                    )

        # trim off last 4 entries; these summarize unaligned features in
        # htseqcounts
        filedata = filedata[0:len(filedata)-5]
        sorted_data = numpy.sort(filedata, order=["featureID"])

        # check if the featurelists are the same:
        try:
            if numpy.all(expressiondata["featureID"] ==
                         sorted_data["featureID"]):
                expressiondata = append_fields(expressiondata,
                                               "%s_rawReads" % filename,
                                               sorted_data["rawReads"],
                                               dtypes="int",
                                               asrecarray=False,
                                               usemask=False)
            else:
                mismatches = numpy.where(expressiondata["featureID"] !=
                                         sorted_filedata["featureID"])[0]
                log.write("ERR: %s contains a different set of feature IDs " \
                          "than in the previously processed files.\n number " \
                          "of IDs in %s: %i\n number of IDs in previous " \
                          "files: %i\n mismatching line IDs: %s\n" \
                          "" % (filename,
                                filename,
                                len(sorted_filedata["featureID"]),
                                len(expressiondata["featureID"]),
                                ", ".join(mismatches)
                                ))
                sys.exit("ERR: %s contains a different set of feature IDs " \
                          "than in the previously processed files.\n number " \
                          "of IDs in %s: %i\n number of IDs in previous " \
                          "files: %i\n mismatching line IDs: %s\n" \
                          "" % (filename,
                                filename,
                                len(sorted_filedata["featureID"]),
                                len(expressiondata["featureID"]),
                                ", ".join(mismatches)
                                ))
            datacolumns.append("%s_rawReads" % filename)

        except TypeError:
            expressiondata = rename_fields(sorted_data,
                                           {"rawReads":"%s_rawReads" \
                                            "" % filename
                                            }
                                           )
        except ValueError:
            pass

    log.write(" %s: Done concatenating data\n" % timeStr())
    log.close()

    if format=="unstruct":
        features = expressiondata["featureID"]
        samples = datacolumns
        # after removing the agnc call for each individual gbg_file, the
        # annotations variable is no longer required. The values in the
        # tuple are used in a hardcoded way later (which I can't be bothtered
        # to rewrite) so I'm including an empty variable here for simplicity.
        # sue me.
        annotations = None
        data = expressiondata[datacolumns]
        data = data.view("int").reshape(data.shape+(-1,))
        return(data, features, annotations, samples)
    else:
        return(expressiondata)

def run_flagstat(bam_file, options):

    ''' Instance local run of samtools flagstat for a bam file

    :param str bam_file: the filename of the bam to summarise

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * the job ID number of the cluster job
        * the name of the file containing the resulting flagstat output

    The *options* instance must contain the keywords:

        * *samtoolspath* - the path to the samtools executable.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.

    The script uses :py:class:`subprocess.Popen` to call an external
    run of samtools flagstat.
    '''

    # open log
    log = open(options.log,"a",0)
    log.write("\tCounting mapped reads in file %s...\n" % bam_file)

    args = [options.samtoolspath, "flagstat", bam_file, ">",
            "%s/%s.flagstatout" % (options.tmpdir, os.path.basename(bam_file))]

    if options.verbose:
        log.write("\t\tcommand line:\n\t\t%s\n" \
              "" % " ".join(args))
        log.write("\t\tJob executing...\n")

    samtools_call = subprocess.Popen(args, stdout=subprocess.PIPE)
    samtools_out = samtools_call.communicate()[0]
    log.write("\t\tcompleted:\n%s\n" % samtools_out)

    log.close()

    return(None, "%s/%s.flagstatout" % (options.tmpdir,
                                        os.path.basename(bam_file)))

def run_flagstat_cluster(c_session, bam_file, options):

    ''' Instance a cluster run of samtools flagstat for a bam file

    :param c_session: a :py:class:drmaa.Session: instance for interacting with
                      the cluster.

    :param str bam_file: the filename of the bam to summarise

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * the job ID number of the cluster job
        * the name of the file containing the resulting flagstat output

    The *options* instance must contain the keywords:

        * *samtoolspath* - the path to the samtools executable.
        * *clustq* - the name of the cluster queue to submit jobs to.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.

    Note that this function includes a SGE-specific job submission parameter
    (in fact, they're quite possibly specific to our individual cluster!).
    Briefly, they are:

        * *-cwd* - specify running the job in the current working directory
        * *-clear* - clear any existing cluster queue settings.
        * *-q 'options.clustq'* - specify the queue to submit the job to.

    '''

    # open log
    log = open(options.log,"a",0)
    log.write("\tCounting mapped reads in file %s...\n" % bam_file)

    c_job = c_session.createJobTemplate()

    c_job.remoteCommand = options.samtoolspath

    args = ["flagstat", bam_file]

    c_job.args = args

    if options.verbose:
        log.write("\t\tdrmaa command line:\n\t\t%s %s\n" \
              "" % (options.samtoolspath," ".join(c_job.args)))

    c_job.outputPath = ":%s" % options.tmpdir
    c_job.errorPath = ":%s" % options.tmpdir

    # pass current working directory (not that this is needed really, but hey!)
    c_job.nativeSpecification = "-cwd"

    # add support for different cluster queue specifications
    c_job.nativeSpecification = "-clear -q '%s' %s" \
                                "" % (options.clustq, c_job.nativeSpecification)


    if options.verbose:
        log.write("\t\tdrmaa output intermediates written to: %s\n" \
                  "" % options.tmpdir)

    c_job.jobEnvironment = os.environ
    jobid = c_session.runJob(c_job)

    log.write("\t\tJob submitted with id: %s\n" % jobid)

    log.close()

    return(jobid, "%s/samtools.o%s" % (options.tmpdir,
                                      jobid))

def deseq_normalization(data, samples):

    ''' Normalize the raw read counts using DESeq normalisation

    :param data: a :py:class:`numpy.ndarray` containing the raw read count
                 summarization data.

    :param list samples: a list of samples (in the order of the columns of the
                         *data* array).

    :return: a tuple containing:

        * a :py:class:`numpy.ndarray` containing the normalized read count
          summarization data.
        * a dictionary keyed on the sample names containing the sample
          normalization factors.

    The form of the normalization is taken from the `DEseq
    <https://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf>`_
    tool implemented in `bioconductor <https://bioconductor.org/>`_. The *data*
    :py:class:`numpy.ndarray` should have columns corresponding to each sample
    in *samples*.

    '''

    # open log
    log = open(options.log,"a",0)
    log.write("\n %s: Running deseq normalization..." % timeStr())

    # calculate per gene geometric mean across replicates
    gm=mst.gmean(data,axis=1)
    gm.mask=numpy.ma.nomask

    # only interested in those with finite geometric mean
    gvs=numpy.where(gm>0)[0]

    # adjust counts by geometric exclude those genes with zero counts
    # This is what DESeq does but is not true representation of the equation
    # it quotes
    dat=numpy.array([(data[i]/gm[i]) for i in gvs if 0 not in data[i]])
    nm=numpy.median(dat,axis=0)

    # calculate normalised count
    norm=numpy.array([(data.transpose()[i]/nm[i]) \
                    for i in numpy.arange(len(nm))],numpy.float).transpose()

    norms={}
    for sample_no in numpy.arange(len(samples)):
        norms[samples[sample_no]]=nm[sample_no]

    return(norm, norms)

def parseFlagstatOut(filename):

    ''' parse the samtools flagstat output for the total number of reads

    :param str filename: the filename of the output from samtools flagstat

    :return: the number of million reads counted by samtools flagstat
    :rtype: float

    '''

    filehandle = open(filename,"r")
    thisdata = filehandle.readlines()
    filehandle.close()


    for line in thisdata:
        try:
            nreads = re.match("^([0-9]+).*mapped.+$", line).group(1)
            break
        except AttributeError:
            pass

    return(float(nreads)/1E6)

def pmrm_normalization(data, all_bams, samples, options):

    ''' Normalize the raw read counts by 'per million mapped reads'

    :param data: a :py:class:`numpy.ndarray` containing the raw read count
                 summarization data.

    :param list all_bams: a list of bam files to extract the total read count
                          from. Should match in length and order the samples
                          listed in *samples*.

    :param list samples: a list of samples (in the order of the columns of the
                         *data* array).

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * a :py:class:`numpy.ndarray` containing the normalized read count
          summarization data.
        * a dictionary keyed on the sample names containing the sample
          normalization factors.

    The *options* instance must contain the keywords:

        * *samtoolspath* - the path to the samtools executable.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.
        * *nocluster* - toggle running the summarization on the cluster, or
          locally.

    Calls either :func:`run_flagstat` or :func:`run_flagstat_cluster` depending
    on whether *options.nocluster* is *True* or *False*. If *options.nocluster*
    is *False* then then *options* must also contain the keyword:

        * *clustq* - the name of the cluster queue to submit jobs to.

    The output from these runs is then parsed with :func:`parseFlagstatOut`.

    '''

    # open log
    log = open(options.log,"a",0)
    log.write("\n %s: Running per-million-reads-mapped normalization..." \
              "" % timeStr())

    # run samtools on each file to get the total number of mapped reads in
    # the ban file. Use the cluster if it is available.
    readcount_files = {}
    if options.nocluster:
        log.write("\n %s: Cluster not enabled: Running samtools flagstat " \
                  "on each bam file in the experiment...\n" \
                  "" % timeStr())

        i=0
        for bam_file in all_bams:
            job, flagstat_outfile = run_flagstat(bam_file, options)
            readcount_files[bam_file] = (flagstat_outfile, i)
            i+=1
    else:

        log.write("\n %s: Cluster enabled: Running samtools flagstat on a " \
                  "cluster node in the %s queue for each bam file in the " \
                  "experiment...\n" % (timeStr(),options.clustq))

        # use drmaa to spread multiple samtools jobs across the cluster.
        # open a drmaa
        c_session = drmaa.Session()
        c_session.initialize()

        # loop throught he files calling a new cluster job for each
        joblist = []
        filelist=[]
        for bam_file in all_bams:
            jobid, flagstat_outfile = run_flagstat_cluster(c_session,
                                                           bam_file,
                                                           options)
            joblist.append(jobid)
            filelist.append(flagstat_outfile)
            readcount_files[bam_file] = (flagstat_outfile, jobid)

        log.write("\n %s: %s jobs launched..." \
                  "" % (timeStr(), os.path.basename(options.samtoolspath)))

        log.close()

        prefix = "%s/%s" % (options.tmpdir,
                            os.path.basename(options.samtoolspath))
        exit_status, sucess, fail, options = monitorDrmaaJobs(c_session,
                                                              joblist, filelist,
                                                              prefix, options,
                                                              delo=False)

        log = open(options.log,"a",0)

        # clearup job metadat and information
        for jobid in joblist:
            c_session.synchronize(jobid,
                                  drmaa.Session.TIMEOUT_WAIT_FOREVER,
                                  True)

        c_session.exit()

    # read files to get readcounts
    norms = {}
    for entry in readcount_files:
        file = readcount_files[entry][0]
        jid = readcount_files[entry][1]
        thissample = None
        sampleid = None
        i=0
        while i<len(samples):
            if re.match(re.sub(".+\/","",file),samples[i]):
                thissample = samples[i]
                sampleid = i
            i+=1

        if options.verbose:
            log.write("\n %s: Reading file %s..." % (timeStr(), file)
                      )
        total_Mreads = parseFlagstatOut(file)
        data[:,sampleid] = data[:,sampleid]/total_Mreads
        norms[file]=total_Mreads

        status = exit_status[jid]
        if status!=drmaa.JobState.UNDETERMINED and status!=drmaa.JobState.FAILED:
            if options.verbose:
                log.write("\n %s: File %s read. Removing file." % (timeStr(), file))
            os.remove(file)

    log.close()

    return(data, norms)

def normalize_data(results, experiment, all_replicates, options):

    ''' Wrapper function to run the various normalization routines

    :param results: a tuple containing:

        * a :py:class:`numpy.ndarray` with the read count summarization
          data.
        * a :py:class:`numpy.ndarray` containing the feature names the
          read count summarization applies to (in the same order as the
          rows of the data array).
        * a list of the samples the read counts were summarized for (in
          the same order as the columns of the data array).

                    These are output as the 0th, 1st and 3rd elements of the
                    tuple from running :func:`read_expression_data` with
                    *format="unstruct"*.

    :param dict experiment: the experiment structure. This is the 1st element
                            of the output tuple from
                            :func:`parse_experiment_data`.

    :param list all_replicates: a list of all the replicates in the experiment
                                (2nd element in the tuple output from
                                :func:`parse_experiment_data`).

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * a :py:class:`numpy.ndarray` containing the normalized read count
          summarization data.
        * a dictionary keyed on the sample names containing the sample
          normalization factors.

    The *options* instance must contain the keyword *norm*. If
    *options.norm='pmrm'* this function calls :func:`pmrm_normalization`, if
    *options.norm='deseq'* this function calls :func:`deseq_normalization`, if
    *options.norm='none'* no normalization is performed, if
    *options.norm='package:<detoolname>'* then no normalization is performed
    here but the DE tool Rscript normalization will be applied further in the
    processing. The *options* instance must contain the keywords appropriate
    for the function being called.

    '''
    # open log
    log = open(options.log,"a",0)

    samples = results[2]
    data = numpy.array(results[0], dtype="float")
    norms = None

    if options.norm == "pmrm":
        data, norms = pmrm_normalization(data, all_replicates, samples, options)
    elif options.norm == "deseq":
        data, norms = deseq_normalization(data, samples)
    elif options.norm == "none":
        log.write("\n %s: performing no normalisation..." % timeStr())
    elif options.norm.split(":")[0] == "package":
        log.write("\n %s: performing package %s normalisation..." % (timeStr(),options.norm.split(":")[-1]))
    else:
        log.write("\n %s: Sorry. '%s' normalisation has yet to be implement. " \
                  "Please try again once you've written it.\n\t\tPerforming no " \
                  "normalisation..." % (timeStr(), options.norm))

    log.close()
    return(data, norms)

def write_R_data_file(results, experiment, options):

    ''' writes text files with data for the tool Rscript routines

    :param results: a tuple containing:

        * a :py:class:`numpy.ndarray` with the read count summarization
          data.
        * a :py:class:`numpy.ndarray` containing the feature names the
          read count summarization applies to (in the same order as the
          rows of the data array).
        * a list of the samples the read counts were summarized for (in
          the same order as the columns of the data array).

                    These are output as the 0th, 1st and 3rd elements of the
                    tuple from running :func:`read_expression_data` with
                    *format="unstruct"*.

    :param dict experiment: the experiment structure. This is the 1st element
                            of the output tuple from
                            :func:`parse_experiment_data`.

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * the filename of the tsv format text data file that the Rscripts
          read as the expression data for differential expression.
        * the filename of the tsv format text data file that the Rscripts
          read as the phenoData for differential expression.

    This data is used by R to contruct expressionSet objects. For details on
    the file format see:

        * http://rss.acs.unt.edu/Rdoc/library/Biobase/html/readExpressionSet.html
        * http://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html

    Basically, we need first need a file with the first column as the feature
    names, then columns for each sample. A header line can contain the column
    names. The input 'results' is a tuple of

       results = (data, features, [samples])

    samples here is optional.

    The *options* instance must contain the keywords:

        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.

    '''

    # open log
    log = open(options.log,"a",0)
    
    R_matrix_filename = "%s/R_exprsFile.txt" % options.tmpdir
    R_pheno_filename = "%s/R_phenoFile.txt" % options.tmpdir

    log.write("\n %s: Writing the matrix file for R (%s, %i columns, %i " \
              "features)..." % (timeStr(),
                                  R_matrix_filename,
                                  len(results[2]),
                                  len(results[1])
                                  )
              )

    exprs_file = open(R_matrix_filename,"w")

    try:
        exprs_file.write("%s\n" % "\t".join(str(val) for val in results[2]))
    except IndexError:
        pass

    i=0
    while i<len(results[1]):
        line_data_string =  "\t".join(str(val) for val in results[0][i])
        exprsFile_line_string = "%s\t%s\n" % (results[1][i], line_data_string)
        exprs_file.write(exprsFile_line_string)
        i+=1

    exprs_file.close()

    log.write("\n %s: Writing the sample information file for R (%s)...\n" \
                  "" % (timeStr(), R_pheno_filename))

    pheno_file = open(R_pheno_filename, "w")

    for sample_name in results[2]:
        for condition in experiment.keys():
            for replicate_file in experiment[condition]:
                if re.match("%s" % os.path.basename(replicate_file),
                            sample_name):
                    pheno_file.write("%s\t%s\n" % (sample_name, condition))

    pheno_file.close()

    log.close()

    return(R_matrix_filename, R_pheno_filename)

def runRDE(rexpfile, rphefile, options, fileout=None, runAGNC=True):

    ''' runs an R-based differential expression analysis for a set of input data

    :param str rexpfile: the filename of the tsv-format text file for the R
                         tool to read. This is the first element of the tuple
                         output from :func:`write_R_data_file`.

    :param str rphefile: the filename of the tsv-format text file for the R
                         tool to read. This is the second element of the tuple
                         output from :func:`write_R_data_file`.

    :param options: a valid :py:class:`optparse.Values` instance.

    :param bool fileout: optional output filename to give to the R tool. If not
                         specified this defaults to the value storred in
                         *options.outfile*.

    :param bool runAGNC: toggle running :func:`run_agnc` on the output data.

    :return: the filename of the output DE calls from the R tool.

    The *options* instance must contain the keywords:

        * *script_file* - the filename of the Rscript to run on the data.
        * *norm* - the normalization toggle for the data.
        * *rpath* - the path to the Rscript executable to use to run
          *script_file*.
        * *outfile* - filename for the Rscript to write its output to.
        * *log* - filename of the the log file to write to.
        * *tmpdir* - the temporary directory to write *restart.pkl* to.
        * *keep_tmpdir* - toggle whether to keep or clean-up the temp dir.
        * *verbose* - boolean toggle for verbose logging.

    The script uses :py:class:`subprocess.Popen` to call an external
    run of `Rscript
    <https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html>`_
    that runs a DE tool on the data. This will get executed ont he cluster.
    Optionally, this function calls :func:`run_agnc` to annotate the returned
    output data.

    '''

    log = open(options.log,"a",0)
    log.write("\n %s: Running R script %s..." % (timeStr(),
                                                 options.script_file))

    if fileout is None:
        fileout = options.outfile

    #fileout_tmp = "%s/%s" % (options.tmpdir, os.path.basename(fileout))
    R_tempfile = "%s/%s_rtemp" % (options.tmpdir, os.path.basename(fileout))
    fileout = "%s/%s_rout" % (options.tmpdir, os.path.basename(fileout))
    
    cmdlist = [options.rpath, options.script_file,rexpfile, rphefile, 
               R_tempfile, options.norm]
    
    if options.verbose:
        log.write("\n\t\commandline: %s" % " ".join(cmdlist))
        log.write("\n\t\tresults stored in temp file %s" % R_tempfile)

    # run R script and capture stderr and stdout. If fails report stderror
    # and exit else only report stdout
    r_call = subprocess.Popen(cmdlist,stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE)

    r_out = r_call.communicate()

    if r_call.returncode>0:
        log.write("\n Failed! stdout and stderr are below:\n %s\n%s" \
                  "" % (r_out[0],r_out[1]))
        raise IOError("R script failed! Can't continue. Halting.")

    log.write("\n\t\t%s" % r_out[0])
    log.write("\n\tFinished\n")

    if runAGNC:
        # add gene descriptions
        log.write("\n %s: Adding gene name and description information to " \
                  "%s..." % (timeStr(), R_tempfile))

        job, R_tempfile = run_agnc(R_tempfile, options)

    # check output file from run_agns() exists and is non-empty
    try:
        sizecheck = os.path.getsize(R_tempfile)
        if sizecheck>0:
            pass
        else:
            # keep all tmpfiles regardless of settings and report error
            options.keep_tmpdir = True
            log.write("\n %s: Failed. Output annotated file is empty. " \
                      "Keeping all temp files, available in %s... " \
                      "" % (timeStr(), options.tmpdir)
                      )
            sys.exit("\n %s: Failed. Output annotated file is empty. " \
                      "Keeping all temp files, available in %s... " \
                      "" % (timeStr(), options.tmpdir)
                      )
    except:
        # keep all tmpfiles regardless of settings and report error
        options.keep_tmpdir = True
        log.write("\n %s: Failed. Output annotated file doesn't exist. " \
                  "Keeping all temp files, available in %s... " \
                  "" % (timeStr(), options.tmpdir)
                  )
        sys.exit("\n %s: Failed. Output annotated file doesn't exist. " \
                  "Keeping all temp files, available in %s... " \
                  "" % (timeStr(), options.tmpdir)
                  )

    log.write("\n %s: writing final output file: %s" % (timeStr(), fileout))

    out_stat = os.stat(R_tempfile)

    filehandle = open(R_tempfile,'r')
    filedata = filehandle.read()
    filehandle.close()

    filehandle = open(fileout,'w')
    filehandle.write('#\n# This file created with command: %s %s %s %s " \
                     "%s' % (options.rpath, options.script_file, rexpfile,
                               rphefile, fileout)
                     )
    filehandle.write('\n# Output created at: %s\n' \
                     "" % timeStr()
                     )
#                     timeAndDateStr()
#                     )

    filehandle.write(filedata)
    filehandle.close()

    log.close()

    return(fileout)

def writeCheckpoint(update, options, checkpointdata=None):

    """ record restart/checkpointing information

    :param dict update: dictionary of run data and results to update/store in a
                        cPickle file.

    :param options: a valid :py:class:`optparse.Values` instance.

    :param dict checkpointdata: a pre-existing dictionary of results to update
                                with the unformation in *update*.

    The *options* instance must contain the keywords:

        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.

    uses :py:class:`cPickle.dump`  to write a file called *restart.pkl* to
    *options.tmpdir*. This is used by the script to restart a failed or halted
    run if the :option:`--restart` option is provided to the script.

    """

    log = open(options.log,"a",0)

    log.write("\n %s: writing checkpointing information..." % timeStr())
    if checkpointdata is None:
        checkpointdata={}

    if type(update) is not dict:
        log.write("Error: checkpoint update data is not in dictionary format!")
        sys.exit()
    else:
        for key in update:
            if key in checkpointdata:
                if options.verbose:
                    log.write("\n\t\tDuplicate checkpoint key detected: %s " \
                              "(overwritting...)" % key)
                checkpointdata[key] = update[key]
            else:
                checkpointdata[key] = update[key]

    file = open("%s/restart.pkl" % options.tmpdir,"w")
    cPickle.dump(checkpointdata, file)
    file.close()

    return(checkpointdata)

def runClusterBootstrap(c_session, i, b, boot_restart_file, options):

    """ Instances a new run of *generic_wrapper.py* script on the cluster

    :param c_session: a :py:class:`drmaa.Session` instance for interacting with
                      the cluster.

    :param int i: index of spawned cluster job cluster. Max value is defined by
                  :option:`--nbootstrap` of the master *generic_wrapper.py*
                  run and is used to define cluster output log filenames.

    :param int b: the number of bootstrap iterations to run in this spawned
                  job. This option is passed as :option:`--nbootstrap`
                  to the spawned slave *generic_wrapper.py* run.

    :param str boot_restart_file: filename of the current restart file.

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * the cluster job ID
        * the filename of the cluster output log
        * the filename of the Rscript output results file
        * the filename of the spawned *generic_wrapper.py* logfile

    This function is called by :func:`clusterBootstrap`. The script uses
    :py:class:`subprocess.Popen` to call an external run of
    *generic_wrapper.py* but with the :option:`--bootstrapslave` set
    to slave it to the master run. The restart file specified
    with *boot_restart_file* is loaded by the spawned job and passes the script
    parameters tot he run. This ensures all bootstrap run use the correct
    parameters and, as a bonus, we don't need to calcualte things twice!

    The *options* instance must contain the keywords:

        * *datapath* - path of the data files to parse.
        * *feature_annot* - the path to the `gff
          <https://www.sanger.ac.uk/resources/software/gff/>`_ file to use for
          gene summarization.
        * *outfile* - the output filename for the spawned cluster run. This is
          passed as :option:`--outfile` to the spawned *generic_wrapper.py* run.
        * *script_file* - the filename of the Rscript to run on the data. This
          is passed as :option:`--scriptfile` to the spawned
          *generic_wrapper.py* run.
        * *clustq* - the name of the cluster queue to submit jobs to.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.

    Note that this function includes a SGE-specific job submission parameter
    (in fact, they're quite possibly specific to our individual cluster!).
    Briefly, they are:

        * *-cwd* - specify running the job in the current working directory
        * *-clear* - clear any existing cluster queue settings.
        * *-q 'options.clustq'* - specify the queue to submit the job to.

    """

    # open log
    log = open(options.log,"a",0)
    log.write("\n\n %s: Cluster bootstrapping enabled. Running job "\
              "%05i\n" % (timeStr(), i+1)
              )

    c_job = c_session.createJobTemplate()

    # run itself!
    #c_job.remoteCommand = "%s %s" % (sys.executable, sys.argv[0])
    c_job.remoteCommand = "%s" % sys.executable

    fileout_details = os.path.splitext(os.path.basename(options.outfile))
    thisjob_fileout = "%s/%s_%04i%s" % (options.tmpdir,
                                        fileout_details[0],
                                        i+1,
                                        fileout_details[1])
    log_details = os.path.splitext(os.path.basename(options.log))
    thisjob_log = "%s/%s_%04i%s" % (options.tmpdir,
                                    fileout_details[0],
                                    i+1,
                                    log_details[1])

    args = [sys.argv[0],
            "-d", options.datapath,
            "-a", options.feature_annot,
            "-o", thisjob_fileout,
            "-l", thisjob_log,
            "-r", options.script_file,
            "-b", str(b),
            "--tmpdir", options.tmpdir,
            "--bootstrapslave"]

    c_job.args = args

    if options.verbose:
        log.write("\t\tdrmaa command line:\n\t\t%s %s\n" \
              "" % (c_job.remoteCommand, " ".join(c_job.args)))

    c_job.outputPath = ":%s" % options.tmpdir
    c_job.errorPath = ":%s" % options.tmpdir

    # pass current working directory (not that this is needed really, but hey!)
    c_job.nativeSpecification = "-cwd"

    # add support for different cluster queue specifications
    c_job.nativeSpecification = "-clear -q '%s' %s" \
                                "" % (options.clustq, c_job.nativeSpecification)

    if options.verbose:
        log.write("\t\tdrmaa output intermediates written to: %s\n" \
                  "" % options.tmpdir)

    c_job.jobEnvironment = os.environ
    jobid = c_session.runJob(c_job)

    log.write("\t\tJob submitted with id: %s\n" % jobid)

    log.close()

    return(jobid, "%s/generic_wrapper.py.o%s" % (options.tmpdir, jobid),
           "%s/%s" % (options.tmpdir, thisjob_fileout), thisjob_log)

def clusterDistribute(nb, nj, options):

    """ calculates optimum bootstrap dist over the specified cluster nodes

    :param int nb: total number of bootstraps

    :param int nj: total number of cluster nodes to use

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a list of tuples each containing:

        * the running total job count
        * the number of of bootstraps to run in that job

    The *options* instance must contain the *log* keywords containing the
    filename of the the log file to write to. This function is called by
    :func:`clusterBootstrap`.

    """

    log = open(options.log,"a",0)
    log.write("\n %s: Calculating the optimum way of spreading %i bootstrap " \
              "runs over %i nodes..." % (timeStr(), nb, nj))

    remaining_nb = float(nb)
    remaining_nj = float(nj)

    bootstrapout=[]
    total_nb = 0
    while remaining_nb>0:

        max_nb_job = math.ceil(remaining_nb/remaining_nj)
        nj_max_nb_job = math.floor(remaining_nb/max_nb_job)
        total_nb = total_nb+(nj_max_nb_job*max_nb_job)

        log.write("\n\t\tAllocating %.0f runs to %.0f nodes (%.0f runs " \
                  "assigned)." % (max_nb_job, nj_max_nb_job, total_nb))

        bootstrapout.append((nj_max_nb_job, max_nb_job))
        remaining_nb = float(remaining_nb-(nj_max_nb_job*max_nb_job))
        remaining_nj = float(remaining_nj-nj_max_nb_job)

    log.close()

    return(bootstrapout)

def clusterBootstrap(options):

    """ setup and run bootstrap iterations on the cluster

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: a tuple containing:

        * a list of cluster logfiles - one for each spawned job
        * a list of output sqlite databases containing DE results - one foe
          each spawned job

    This scrips takes in input options for the bootstrap runs and uses them to
    spawn a series of cluster runs of *generic_wrapper.py* (in slave mode) to
    actually do the DE calculations.

    The *options* instance must contain the keywords:

        * *nbootstrap* - number of bootstraps to run.
        * *njobs* - the max number of cluster jobs to spawn.
        * *restart* - toggle whether to use the *restart.pkl* file to start a
          run - this is used by the runs spawned by this function.
        * *nocluster* - toggle runnign the bootstraps on the cluster.
        * *bootmaster* - internally set variable denoting the input as the
          master script.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to write *restart.pkl* to.
        * *verbose* - boolean toggle for verbose logging.

    The *options* instance must also contain all the keywords required by
    :func:`runClusterBootstrap`, :func:`clusterDistribute` and
    :func:`monitorDrmaaJobs`, all of which are called by this function and are
    passed the input :py:class:`optparse.Values` instance.

    """

    log = open(options.log,"a",0)
    log.write("\n %s: running in cluster-bootstrap mode with %i " \
              "iterations" % (timeStr(), int(options.nbootstrap)))
    log.write("\n %s: writing bootstrap data pickle..." % timeStr())

    # read existing restart information
    file = open("%s/restart.pkl" % options.tmpdir,"r")
    bootstrap_restardic = cPickle.load(file)
    file.close()

    # modify bootstrap options; set bootstrap slave, --nocluster
    # and get the load distributions
    bootstrap_restardic["options"].bootmaster = False
    bootstrap_restardic["options"].nocluster = True
    bootstrap_restardic["options"].restart = True

    # write new restart pkl
    boot_restart_file = "%s/bootruns_restart.pkl" % options.tmpdir
    file = open(boot_restart_file,"w")
    cPickle.dump(bootstrap_restardic, file)
    file.close()

    # use drmaa to spread multiple jobs across the cluster.
    c_session = drmaa.Session()
    c_session.initialize()

    jobloads = clusterDistribute(int(options.nbootstrap),
                                int(options.njobs),
                                options)

    # loop through the jobloads:
    i=0
    joblist = []
    boostrap_files = []
    boostrap_logs = []
    boostrap_dbs = []
    for jobload in jobloads:
        j=0
        while j<jobload[0]:
            results = runClusterBootstrap(c_session, i, int(jobload[1]),
                                          boot_restart_file, options)

            joblist.append(results[0])
            boostrap_files.append(results[1])
            boostrap_logs.append(results[2])
            boostrap_dbs.append(results[3])
            i+=1
            j+=1

    log.write("\n %s: %s jobs launched..." \
              "" % (timeStr(), os.path.basename(sys.argv[0]))
              )

    prefix = "%s/%s" % (options.tmpdir, os.path.basename(sys.argv[0]))
    exit_status, sucess, fail, options = monitorDrmaaJobs(c_session, joblist,
                                                          boostrap_dbs, prefix,
                                                          options, delo=False,
                                                          dele=False)

    # clearup job metadata and information
    for jobid in joblist:
        c_session.synchronize(jobid,
                              drmaa.Session.TIMEOUT_WAIT_FOREVER,
                              True)

    c_session.exit()

    log.close()
    return(boostrap_logs, boostrap_dbs)

def results2sqlite(features, results_files, options):

    """ converts the text output from n bootstrap runs into one sqlite db

    :param list features: a list of feature identified the DE has been
                          calculated for.

    :param list results_files: a list of DE output text files with a common
                              format. Output from :func:`runRDE`.

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: the filename of the resulting sqlite database, containing the
             results from all the files listed in results_files.

    The *options* instance must contain the keywords:

        * *log* - filename of the the log file to write to.
        * *verbose* - boolean toggle for verbose logging.

    This function takes the results from multiple DE textfiles and concatenates
    them into a single sqlite database that is written to the same directory as
    the input file. In this case, this should be the path given by
    *options.tmpdir* passed to :func:`runRDE`. The database structure is
    described :ref:`here <genwrapout>`.

    """

    # open log
    log = open(options.log,"a",0)

    dbfile = "%s/%s" % (os.path.dirname(results_files[0]), 
                        re.sub("_rout","",os.path.basename(results_files[0])))

    if os.path.exists(dbfile):
            log.write("\tremoving pre-existing version of %s... \n" \
                      "" % dbfile
                      )
            os.remove(dbfile)

    log.write("\n %s: merging bootstrap results to sqlite database (%s)" \
              "" % (timeStr(), dbfile))

    # make an sqlite db to hold the results (with foreign key support)
    if options.verbose:
        log.write("\n\t\tbuilding schema...")

    conn = sqlite3.connect(dbfile)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys=ON")

    # create features table
    c.execute("CREATE TABLE features(" \
                "id INTEGER PRIMARY KEY, "\
                "featureID TEXT, " \
                "name TEXT, " \
                "desc BLOB)"
              )

    # create bootstrap iteration table
    c.execute("CREATE TABLE bootstraps(" \
                "id INTEGER PRIMARY KEY, " \
                "comments BLOB, " \
                "deruntime_sec REAL)"
              )

    # define cols for foreign key links
    col_str = "id INTEGER PRIMARY KEY, " \
              "featureid INTEGER, " \
              "bsid INTEGER, "

    # define foreign key links
    fkey_str = "FOREIGN KEY(featureid) REFERENCES features(id), " \
               "FOREIGN KEY(bsid) REFERENCES bootstraps(id)"

    # get column headers
    cols = None

    file = open(results_files[0],"r")
    for line in file:
        if re.match("^#",line):
            pass
        else:
            cols = re.split("\t", re.sub("\.","",re.sub('"',"",line.strip())))
            break

    file.close()
        
    # construct create table columns
    col_str = "%s%s REAL, %s" % (col_str," REAL, ".join(cols[1:]), fkey_str)
    
    # create DEresutls table - shouldn't realy use %s here, but screw it!
    c.execute("CREATE TABLE DEresults(%s)" % col_str)

    # first insert features
    if options.verbose:
        log.write("\n\t\tloading %i features..." % len(features))

    # features needs to be converted to a list of tuples first
    feature_tups = list((l,) for l in features)
    c.executemany('INSERT INTO features VALUES (NULL,?,"","")', feature_tups)

    # get features primary-key dictionary
    c.execute("SELECT featureid, id FROM features")
    feature_dic = dict(c.fetchall())#

    if options.verbose:
        log.write("\n\t\tmunging results from files...")

    i=0
    filedata=[]
    while i<len(results_files):
        if options.verbose:
            filehandle = open(results_files[i],"r")
            filelines = filehandle.readlines()
            filehandle.close()
            log.write("\n\t\t\t%s (%i lines)" % (results_files[i],
                                                 len(filelines))
                      )

        bsid = None
        file = open(results_files[i],"r")
        j=0
        ncom=0
        ndat=0
        comments=""
        for line in file:
            if re.match("^#",line):
                comments = "%s%s" % (comments, line)
                ncom+=1
            elif j==ncom:
                db_tup = (comments, 0.0)
                c.execute("INSERT INTO bootstraps VALUES(NULL, ?, ?)",db_tup)
                bsid = c.lastrowid
            else:
                linedata = re.split("\t", re.sub('"',"",line.strip()))

                k=1
                datalist=[feature_dic[linedata[0]], bsid]
                while k<len(linedata):
                    datalist.append(linedata[k])
                    k+=1

                filedata.append(tuple(datalist))
                ndat+=1
            j+=1
        file.close()
        i+=1

    if options.verbose:
        log.write("\n\t\tinserting data...")

    de_insertstr = "INSERT INTO DEresults VALUES (NULL,%s)" \
                   "" % re.sub(",$","","?,"*(len(cols)+1))
    
    print de_insertstr
    c.executemany(de_insertstr, filedata)

    conn.commit()
    conn.close()

    return(dbfile)

def concatSqliteFiles(options):

    ''' Concatenates many sqlite files with the same schema

    :param options: a valid :py:class:`optparse.Values` instance.

    :return: the output filename for the concatenated sqlite database.

    The *options* instance must contain the keywords:

        * *outfile* - the output filename for the final file.
        * *log* - filename of the the log file to write to.
        * *tmpdir*- the temporary directory to search for *.db* files in.
        * *verbose* - boolean toggle for verbose logging.

    This function takes in multiple sqlite database files with the same
    structure and concatenates them into a single output database. It assumes
    that the intermediaty database files are located in the path listed in
    *options.tmpdir* and also assumes they have the unique extension *.db*.

    '''

    # open log
    log = open(options.log,"a",0)

    alltmpfiles = os.listdir(options.tmpdir)
    dbfiles=[]
    for file in alltmpfiles:
        if os.path.splitext(file)[1]==".db":
            dbfiles.append(file)
    ## TODO - check that dbfiles actually found some files. Otherwise crashes
    ##        later...
    ##      - in fact this should check that the number of dbfiles matches the
    ##        number of reps

    dbfile = options.outfile
    log.write("\n %s: merging bootstrap results to final sqlite database (%s)" \
              "" % (timeStr(), dbfile))

    if os.path.exists(dbfile):
        log.write("\n %s: sqlite database file (%s) already exists! " \
                  "Overwriting... " % (timeStr(), dbfile))
        os.remove(dbfile)

    if options.verbose:
        log.write("\n\t\tbuilding schema...")

    # make an sqlite db to hold the results (with foreign key support)
    conn = sqlite3.connect(dbfile)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys=ON")

    # create features table
    c.execute("CREATE TABLE features(" \
                "id INTEGER PRIMARY KEY, "\
                "featureID TEXT, " \
                "name TEXT, " \
                "desc BLOB)"
              )

    # create logfile table
    c.execute("CREATE TABLE bslogs(" \
                "id INTEGER PRIMARY KEY, " \
                "log BLOB)"
              )

    # create bootstrap iteration table
    c.execute("CREATE TABLE bootstraps(" \
                "id INTEGER PRIMARY KEY, " \
                "bsinlogid INTEGER, " \
                "comments BLOB, " \
                "deruntime_sec REAL, " \
                "logid INTEGER, " \
                "FOREIGN KEY(logid) REFERENCES bslogs(id))"
              )

    # define cols for foreign key links
    col_str = "id INTEGER PRIMARY KEY, " \
              "featureid INTEGER, " \
              "bsid INTEGER, "

    # define foreign key links
    fkey_str = "FOREIGN KEY(featureid) REFERENCES features(id), " \
               "FOREIGN KEY(bsid) REFERENCES bootstraps(id)" \
               ""

    # connect to the first file
    c.execute("ATTACH '%s/%s' AS db2" % (options.tmpdir, dbfiles[0]))

    # get columns
    c.execute("PRAGMA db2.table_info(DEresults)")
    coldata = c.fetchall()
    colids = []
    i=0
    while i<len(coldata):
        colids.append(coldata[i][1])
        i+=1

    col_str = "%s%s REAL, %s" % (col_str," REAL, ".join(colids[3:]), fkey_str)

    # create DEresutls table - shouldn't realy use %s here, but screw it!
    c.execute("CREATE TABLE DEresults(%s)" % col_str)

    if options.verbose:
        log.write("\n\t\tadding results and log information from "
                  "itermiadary database files...")

    gfeatures=None
    for file in dbfiles:
        file = "%s/%s" % (options.tmpdir, file)
        if options.verbose:
            log.write("\n\t\t\t%s" % file)

        # read the relevant logfile and put log into the database
        logfile = "%s.log" % os.path.splitext(file)[0]
        filehandle = open(logfile,"r")
        logdata = filehandle.read()
        filehandle.close()
        dbtup = (logdata,)
        c.execute("INSERT INTO bslogs VALUES(NULL,?)", dbtup)
        logid = c.lastrowid

        # connect to the relevant database
        conn2 = sqlite3.connect(file)
        c2 = conn2.cursor()
        c2.execute("PRAGMA foreign_keys=ON")

        # if its the first run, add feature information to the database
        if gfeatures is None:
            c2.execute("SELECT featureID, name, desc FROM features")
            features = c2.fetchall()
            c.executemany("INSERT INTO features VALUES(NULL,?,?,? )", features)
            c.execute("SELECT * FROM features")
            gfeatures = c.fetchall()

        # check that the features all appear int he same order.
        c2.execute("SELECT * FROM features")
        features = c2.fetchall()
        if features!=gfeatures:
            print "uh oh!"
            filehandle=open("featdump.pkl","w")
            cPickle.dump((features, gfeatures),filehandle)
            filehandle.close()

        # then get bs metadata and add the data that corresponds to that
        # bootstrap
        c2.execute("SELECT * FROM bootstraps")
        bsinfo = c2.fetchall()

        for val in bsinfo:
            # first insert bootstrap info
            c.execute("INSERT INTO bootstraps VALUES(NULL,?,?,?,%s)" \
                      "" % logid, val)
            bsid = c.lastrowid

            # get data corresponding to this bootstrap run
            c2.execute("SELECT featureid, %s FROM DEresults where bsid=?" \
                       "" % ", ".join(colids[3:]),
                       (val[0],)
                       )
            c2data = c2.fetchall()

            insert_str = "INSERT INTO DEresults VALUES(NULL,?,%i,%s)" \
                         "" % (bsid, re.sub(",$","","?,"*(len(c2data[0])-1)))
            c.executemany(insert_str, c2data)

    conn.commit()
    conn.close()

    return(dbfile)

if __name__ == '__main__':

    # parse command line options - note that because we're using python 2.6.4,
    # and hence optparse rather than argparse, I'm forcing everything to be an
    # 'option' even if it is required.

    # define options
    optslist = []
    optslist.append({
                     "short":   "-d",
                     "long":    "--datapath",
                     "dest":    "datapath",
                     "action":  "store",
                     "help":    "Path to the data directory of the " \
                                "experiment. Each sub-directory in this path" \
                                "will be treated as a separate condition in " \
                                "the experiment (with the name of the " \
                                "condition matching the name of the " \
                                "directory), and each .bam file (or .gbgout " \
                                "files if the --precounts flag is set) in " \
                                "each directory is a replicate for that " \
                                "condition (BAM file indexes for each file "\
                                "should have the same filename, with .bai " \
                                "post-pended).",
                     "group":    "required",
                     "default":  None,
                     "opt_type": "input_path"
                    })
    optslist.append({
                     "short":    "-a",
                     "long":     "--annotation",
                     "dest":     "feature_annot",
                     "action":   "store",
                     "help":     "Path to the .gff feature annotation file " \
                                 "for the data. This file should match the " \
                                 "feature you are counting the RNA-Seq " \
                                 "expression for.",
                     "group":    "required",
                     "default":  None,
                     "opt_type": "input_path"
                    })
    optslist.append({
                     "short":    "-r",
                     "long":     "--scriptfile",
                     "dest":     "script_file",
                     "action":   "store",
                     "help":     "Name and path of script file to execute " \
                                 "in order to perform differential gene " \
                                 "expression analysis.",
                     "group":    "required",
                     "default":  None,
                     "opt_type": "input_file"
                     })
    optslist.append({
                     "short":    "-o",
                     "long":     "--outfile",
                     "dest":     "outfile",
                     "action":   "store",
                     "help":     "The name (inc. path if different from the " \
                                 "current working directory) of the output " \
                                 "file containing the DE results.",
                     "group":    "required",
                     "default":  None,
                     "opt_type": "output_file"
                    })
    optslist.append({
                     "short":    "-l",
                     "long":     "--logfile",
                     "dest":     "log",
                     "action":   "store",
                     "help":     "The name (inc. path if different from the " \
                                 "current working directory) of the log " \
                                 "file from running this script.",
                     "group":    "required",
                     "default":  None,
                     "opt_type": "output_file"
                    })
    optslist.append({
                     "short":   "-s",
                     "long":    "--species",
                     "dest":    "species",
                     "action":  "store",
                     "help":    "The species used to generate the RNA-Seq " \
                                "data. This is used by the ensembl API to " \
                                "annotate the features found, so it must be " \
                                "recognizable to ensembl. For convenience " \
                                "the default is 'Saccharomyces cerevisiae'.",
                     "group":   "None",
                     "default": 's_cerevisiae'
                    })
    optslist.append({
                     "short":   "-v",
                     "long":    "--verbose",
                     "dest":    "verbose",
                     "action":  "store_true",
                     "help":    "Turn on verbose logging.",
                     "group":   None,
                     "default": False
                    })
    optslist.append({
                "short":    None,
                "long":     "--precounts",
                "dest":     "precounts",
                "action":   "store_true",
                "help":     "If this option is set, the script looks for pre-" \
                            "generated gene count files in  the --datapath. " \
                            "If found, these gene count files are used in " \
                            "place of bam files and gene count summarization " \
                            "is skipped. If this option is set and no gene " \
                            "count files are found then the script will " \
                            "default tot he bam files and will do gene count " \
                            "summarization. If the --savecounts flag is set, " \
                            "with this option and no gene count files are " \
                            "found, summarization will be done and the files " \
                            "saved to --datapath. Defaults to False.",
                "default":  False,
                "group":    "counting"
                })
    optslist.append({
                "short":    None,
                "long":     "--precount_fileext",
                "dest":     "precount_fileext",
                "action":   "store",
                "help":     "The file extension of pre-generated gene count " \
                            "data. Defaults to 'gbgout'.",
                "default":  'gbgout',
                "group":    "counting"
                })
    optslist.append({
                "short":    None,
                "long":     "--savecounts",
                "dest":     "savecounts",
                "action":   "store_true",
                "help":     "If this flag is set, the gene count files from " \
                            "running the count summarization part of the " \
                            "code  will be saved to --datapath",
                "default":  False,
                "group":    "counting"
                })
    optslist.append({
                     "short":        "-f",
                     "long":         "--feature",
                     "dest":         "feature",
                     "action":       "store",
                     "help":         "Sets the level at which you count the " \
                                     " RNA-Seq signal. Options are: gene " \
                                     "(default) & transcript.",
                     "group":        "counting",
                     "default":      "gene_id",
                     "allowed_vals": ["gene_id","transcript_id"]
                    })
    optslist.append({
                     "short":        "-c",
                     "long":         "--count-method",
                     "dest":         "cmeth",
                     "action":       "store",
                     "help":         "Choose the method used for counting " \
                                     " the reads for each feature. Options " \
                                     "are: union (default), istrict, " \
                                     "noempty & cufflinks (not yet used)",
                     "group":        "counting",
                     "default":      "union",
                     "allowed_vals": ["union", "istrict",
                                      "inoempty"]
                    })
    optslist.append({
                     "short":        "-n",
                     "long":         "--norm",
                     "dest":         "norm",
                     "action":       "store",
                     "help":         "Normalization method for the tool " \
                                     "being wrapped. Options are: " \
                                     "package:normtype (uses a normtype " \
                                     "method available in the package), "\
                                     "none (default), pmrm (norm by Per " \
                                     "Million Reads Mapped), fpkm (norm by " \
                                     "Fragments Per Kilobase of exon per " \
                                     "Million fragments mapped), rpkm norm " \
                                     "by Reads Per Kilobase of exon per " \
                                     "Million reads mapped), tmm (norm by " \
                                     "Trimmed Mean of M component) & quartile.",
                     "group":        "normalization",
                     "default":      "package:normtype"
                    })
    optslist.append({
                     "short":    None,
                     "long":     "--gbgfile",
                     "dest":     "gbgpath",
                     "action":   "store",
                     "help":     "The full path to the group_by_gene.pl " \
                                 "script that this script will call to " \
                                 "run htseq-count. The default location " \
                                 "for group_by_gene.pl is in the same " \
                                 "directory as this script.",
                     "group":    None,
                     "default":  "Bootstrapping/group_by_gene.pl",
                     "opt_type": "input_path"
                    })
    optslist.append({
                     "short":    None,
                     "long":     "--agncfile",
                     "dest":     "agncpath",
                     "action":   "store",
                     "help":     "The full path to the " \
                                 "add_gene_name_column.pl script that this " \
                                 "script will call to annotate the output " \
                                 "from group_by_gene.pl. The default " \
                                 "location for add_gene_name_column.pl is " \
                                 "in the same directory as this script.",
                     "group":    None,
                     "default":  "Bootstrapping/add_gene_name_column.pl",
                     "opt_type": "input_path"
                    })
    optslist.append({
                     "short":    None,
                     "long":     "--Rpath",
                     "dest":     "rpath",
                     "action":   "store",
                     "help":     "Specify path to 'Rscript' executable. " \
                                 "Defaults to /sw/opt/R/3.2.2/bin/Rscript",
                     "group":    None,
                     "default":  "/usr/bin/Rscript",
                     "opt_type": "input_file"
                     })
    optslist.append({
                     "short":    None,
                     "long":     "--perlpath",
                     "dest":     "perlpath",
                     "action":   "store",
                     "help":     "Specify path to perl executable." \
                                 "Defaults too /usr/bin/perl",
                     "group":    None,
                     "default":  "/usr/bin/perl",
                     "opt_type": "input_file"
                     })
    optslist.append({
                     "short":   None,
                     "long":    "--tmpdir",
                     "dest":    "tmpdir",
                     "action":  "store",
                     "help":    "The full path of the temp dir you with to " \
                                "use for storing intermediate files. The " \
                                "default location for this is a dir called " \
                                "'.DEtemp' in the current directory.",
                     "group":   None,
                     "default": "%s/.DEtemp" % os.getcwd(),
                    })
    optslist.append({
                     "short":   None,
                     "long":    "--keep-tmpdir",
                     "dest":    "keep_tmpdir",
                     "action":  "store_true",
                     "help":    "Toggle whether to keep tmpdir on " \
                                "successful completion. Default is to not " \
                                "keep tmpdir.",
                     "group":   None,
                     "default": False
                     })
    optslist.append({
                     "short":   None,
                     "long":    "--nocluster",
                     "dest":    "nocluster",
                     "action":  "store_true",
                     "help":    "Enable cluster support through " \
                                "drmaa-python. The default is True",
                     "group":   None,
                     "default": False,
                    })
    optslist.append({
                     "short":   None,
                     "long":    "--clustq",
                     "dest":    "clustq",
                     "action":  "store",
                     "help":    "Set the q to use to run cluster jobs. " \
                                "Default is 64bit-pri.q.",
                     "group":   None,
                     "default": "64bit-pri.q",
                    })
    optslist.append({
                     "short":    None,
                     "long":     "--samtoolspath",
                     "dest":     "samtoolspath",
                     "action":   "store",
                     "help":     "The full path to the samtools instalation " \
                                 "you want to use to manipulate sam/bam " \
                                 "files. The default location for samtools " \
                                 "is /local/bin/samtools",
                     "group":    None,
                     "default":  "/usr/bin/samtools",
                     "opt_type": "input_path"
                    })
    optslist.append({
                     "short":   None,
                     "long":    "--restart",
                     "dest":    "restart",
                     "action":  "store_true",
                     "help":    "Restart a failed job. Will attempt to run " \
                                "process after the last successfully " \
                                "completed step.",
                     "group":   None,
                     "default": False
                    })
    optslist.append({
                     "short":    "-i",
                     "long":     "--includelist",
                     "dest":     "inc_reps",
                     "action":   "store",
                     "help":     "Filename with the filenames of specific " \
                                 "replicates that must be included in the data",
                     "group":    "selection",
                     "default":  None,
                     "opt_type": "input_file"
                     })
    optslist.append({
                    "short":    "-e",
                    "long":     "--excludelist",
                    "dest":     "exc_reps",
                    "action":   "store",
                    "help":     "Filename with the filenames of specific" \
                                "replicates that must be excluded from the data",
                     "group":    "selection",
                    "default":  None,
                    "opt_type": "input_file"
                    })
    optslist.append({
                "short":    "-k",
                "long":     "--kreps",
                "dest":     "kreps",
                "action":   "store",
                "help":     "The number of replicates to use. Replicates " \
                            "listed in --excludelist are subtracted first. " \
                            "Then, if this is None (default) or the number " \
                            "is greater than remaining number of replicates,"\
                            " then all the replicates are used. If this " \
                            "value is less than the remaining number of " \
                            "replicates, then this number of replicates is " \
                            "selected at random for each run.",
                "default":  None,
                "group":    "selection"
                })
    optslist.append({
                "short":    "-b",
                "long":     "--nbootstrap",
                "dest":     "nbootstrap",
                "action":   "store",
                "help":     "The number of bootstrap iterations to run. " \
                            "Default is None",
                "default":  None,
                "group":    "selection"
                })
    optslist.append({
                "short":    None,
                "long":     "--bootstrapmaster",
                "dest":     "bootmaster",
                "action":   "store_true",
                "help":     "Is this a master script initiating a bootstrap? " \
                            "Default is False. This parameter should never " \
                            "be changed by the user, it is set internally.",
                "default":  True,
                "group":    "selection"
                })
    optslist.append({
                "short":    None,
                "long":     "--bootstrapslave",
                "dest":     "bootslave",
                "action":   "store_true",
                "help":     "Is this a slave script for a cluster bootstrap? " \
                            "Default is False. This parameter should never " \
                            "be changed by the user, it is set internally.",
                "default":  False,
                "group":    "selection"
                })
    optslist.append({
                "short":    "-j",
                "long":     "--njobs",
                "dest":     "njobs",
                "action":   "store",
                "help":     "The number of cluster jobs to spawn if not " \
                            "--nocluster. Default is 400",
                "default":  400,
                "group":    "selection"
                })

    # usage string
    usage = """
%prog -d|--datapath <path> -a|--annotation <file> -r|--scriptfile <file>
-o|--outfile <file> -l|--logfile <file> [-s|species <str>] [--precounts]
[--precount_fileext <str>] [--savecounts] [-f|--feature <str>]
[-c|--count-method <str>] [-n|--norm <str>] [-k|--kreps <int>]
[-b|--nbootstrap <int>] [--includelist <file>] [--excludelist <file>]
[--gbgfile <path><filename>] [--agncfile <path><filename>][--Rpath <path>]
[--perlpath <path>] [--tmpdir <path>] [--keep-tmpdir] [--samtoolspath <path>]
[--njobs <int>] [--clustq <str>] [--nocluster] [--tmpdir <path>]
[--keep-tempdir] [--restart] [--version] [-v|--verbose] [--help]
"""

    parser = OptionParser(usage, version="%prog v"+version_string)

    #define option groups
    reqgroup = OptionGroup(parser, "Required Arguments")
    csgroup = OptionGroup(parser,
                          "Feature Counting Options",
                          "This script uses htseq-count, implementd by the " \
                          "script 'group_by_gene.pl' to summarize RNA-Seq " \
                          "read counts for the specified genomic features. " \
                          "See group_by_gene.pl --man for more details on " \
                          "this script. See http://tinyurl.com/cg88yso for " \
                          "more on htseq-count. Note: group_by_gene.pl " \
                          "doesn't currently support choosing different " \
                          "feature types, or different counting methods, so " \
                          "these options do nothing at the moment!")
    normgroup = OptionGroup(parser, "Normalization Options",
                            "Many Differential Expression tools are " \
                            "designed to use a specific normalization. This " \
                            "script keeps this option, but where possible " \
                            "allows 'standard' normalizations so that the " \
                            "impact of changing the normalization can be " \
                            "examind, and so that the results from different " \
                            "DE tools can be compared on an even footing.")
    selgroup = OptionGroup(parser, "Selection and Bootstrap Options",
                            "These options concern subselecting the total " \
                            "number of replicates and running bootstrap " \
                            "iterations of tools. They allow you to specify " \
                            "individual replicates you want to include, or " \
                            "exclude, as well as how many random replicates " \
                            "to select from the remaining set, and how many " \
                            "bootstrap iterations to perform.")

    # add options to the groups
    for val in optslist:
        try:
            if val["group"]=="required":
                reqgroup.add_option(val["short"], val["long"],
                                    action=val["action"], dest=val["dest"],
                                    help=val["help"], default=val["default"])
            elif val["group"]=="counting":
                if val["short"] is not None:
                    csgroup.add_option(val["short"], val["long"],
                                       action=val["action"],
                                       dest=val["dest"], help=val["help"],
                                       default=val["default"])
                else:
                    csgroup.add_option(val["long"], action=val["action"],
                                       dest=val["dest"], help=val["help"],
                                       default=val["default"])

            elif val["group"]=="normalization":
                normgroup.add_option(val["short"], val["long"],
                                     action=val["action"], dest=val["dest"],
                                     help=val["help"], default=val["default"])
            elif val["group"]=="selection":
                if val["short"] is not None:
                    selgroup.add_option(val["short"], val["long"],
                                        action=val["action"], dest=val["dest"],
                                        help=val["help"],
                                        default=val["default"])
                else:
                    selgroup.add_option(val["long"], action=val["action"],
                                        dest=val["dest"], help=val["help"],
                                        default=val["default"])
            else:
                if val["short"] is not None:
                    parser.add_option(val["short"], val["long"],
                                      action=val["action"], dest=val["dest"],
                                      help=val["help"], default=val["default"])
                else:
                    parser.add_option(val["long"], action=val["action"],
                                      dest=val["dest"], help=val["help"],
                                      default=val["default"])
        except TypeError, e:
            print val
            raise(e)

    # add groups to the parser
    parser.add_option_group(reqgroup)
    parser.add_option_group(csgroup)
    parser.add_option_group(normgroup)
    parser.add_option_group(selgroup)

    # parse arguments
    options, arguments = parser.parse_args()

    # ==========
    # Start with housekeeping, option parsing
    # ==========
    restart_dic = None

    if options.bootslave:

        # load generic data if running as a bootstrap slave
        file = open("%s/bootruns_restart.pkl" % options.tmpdir, "r")
        restart_dic = cPickle.load(file)
        file.close()

        # store run specfic options
        thisjob_fileout = options.outfile
        thisjob_b = options.nbootstrap
        thisjob_log = options.log

        # overwrite options
        options = restart_dic["options"]

        options.masterworkingdir = os.path.abspath(options.tmpdir)
        options.masterfileout = os.path.abspath(thisjob_fileout)
        
        if options.verbose:
            options.tmpdir = os.path.abspath(tempfile.mkdtemp(dir=options.tmpdir))
        else:
            # overwrite tmpdir if there is an appropriate local env var
            if "SGE_TEMP" in os.environ.keys():
                options.tmpdir = "%s/%s" % (os.environ["SGE_TEMP"],
                                            os.path.basename(options.tmpdir))
            elif "TMPDIR" in os.environ.keys():
                options.tmpdir = "%s/%s" % (os.environ["TMPDIR"],
                                            os.path.basename(options.tmpdir))

        # restore run-specific options
        options.outfile = "%s/%s" % (os.path.dirname(options.tmpdir),
                                     os.path.basename(thisjob_fileout))

        options.nbootstrap = thisjob_b
        options.bootmaster=False
        options.bootslave = True
        options.log = os.path.abspath(thisjob_log)

        # open the logfile and write the basic header information
        write_log_header(options, optslist, version_string=version_string,
                         imported_modules=imported_modules)

        log = open(options.log,"a",0)
        log.write("\n\n %s: Running as a cluster bootstrap iteration. " \
                  "Loading pre-processed details..." % timeStr())

        createNewTempdir(options, remove=True, create=True)
        os.chdir(options.tmpdir)

        if options.verbose:
            #log.write("\t\tworking in: %s" % os.environ["TMPDIR"])
            log.write("\n\t\ttmpdir is: %s" % options.tmpdir)

        log.close()

    elif options.restart:
        file = open("%s/restart.pkl" % options.tmpdir,"r")
        restart_dic = cPickle.load(file)
        file.close()
        options = restart_dic["options"]

        log = open(options.log,"a",0)
        log.write("\n\n %s: Restarting run from file at restart position " \
                  "%i..." % (timeStr(), restart_dic["restarpos"]))
        log.close()

    else:
        # check that all required options have been provided
        parse_required_options(options, optslist)

        # check all values are allowed, or default to, well, the defaults!
        parse_allowed_values(options, optslist)

        # check that all options are of the right type (if specified)
        options = parse_option_type(options, optslist)

        # open the logfile and write the basic header information
        write_log_header(options, optslist, version_string=version_string)

        # remove/create tmpdir
        createNewTempdir(options, remove=True, create=True)

    # ==========
    # Parse the experiment directory structure
    #     - allow for restart
    # ==========
    experiment = None
    all_replicates = None
    includelist = None

    if restart_dic is None:
        experiment_data = parse_experiment_data(options)
        restart_update = {"restarpos":1,
                          "options":options,
                          "experiment":experiment_data[0],
                          "all_replicates":experiment_data[1],
                          "includelist":experiment_data[2],
                          "do_genecounts":experiment_data[3]}
        restart_dic = writeCheckpoint(restart_update, options)
    else:
        log = open(options.log,"a",0)
        log.write("\n %s: loaded experiment details" % timeStr())
        log.close()

    experiment = restart_dic["experiment"]
    all_replicates = restart_dic["all_replicates"]
    includelist = restart_dic["includelist"]
    do_genecounts = restart_dic["do_genecounts"]

    # ==========
    # Map the reads in bam files to Genes using group_by_gene script
    #     - allow for restart
    # ==========
    gene_count_files = None
    if restart_dic["restarpos"]<2:
        if do_genecounts:
            gene_count_files = get_gene_counts(all_replicates, options)
        else:
            gene_count_files = all_replicates

        restart_update={"gene_count_files":gene_count_files,
                        "restarpos":2}

        restart_dic = writeCheckpoint(restart_update, options,
                                      checkpointdata=restart_dic)
    else:
        gene_count_files = restart_dic["gene_count_files"]

        log = open(options.log,"a",0)
        log.write("\n %s: loaded gene_count file details" % timeStr())
        log.close()

    # ==========
    # Read the gene count files and aggregate them to a data structure
    #     - allow for restart
    # ==========
    exprs_data = None
    if restart_dic["restarpos"]<3:
        # read datafiles
        exprs_data = read_expression_data(gene_count_files, options)

        restart_update = {"exprs_data":exprs_data,
                          "restarpos":3}

        restart_dic = writeCheckpoint(restart_update, options,
                                      checkpointdata=restart_dic)
    else:
        exprs_data = restart_dic["exprs_data"]

        log = open(options.log,"a",0)
        log.write("\n %s: loaded gene expression data" % timeStr())
        log.close()

    # ==========
    # Branch the code for performing a single calculation, or performing a
    # replicate bootstrap.
    # ==========

    if options.nbootstrap is None:
        options.bootmaster = False
        options.bootslave = False
    elif options.nocluster:
        options.bootmaster = False
        options.bootslave = True
    elif not options.nocluster:

        # construct and launch the cluster bootstrap jobs...
        outfile = None
        if restart_dic["restarpos"]<8:
            outfiles = clusterBootstrap(options)
            restart_update = {"outfiles":outfiles,
                              "restarpos":8}

            restart_dic = writeCheckpoint(restart_update, options,
                              checkpointdata=restart_dic)
        else:
            outfiles = restart_dic["outfiles"]

            log = open(options.log,"a",0)
            log.write("\n %s: loaded intermediary databases" % timeStr())
            log.close()

        # concatenate the job databases to a single sqlite database
        if restart_dic["restarpos"]<9:
            masterout = concatSqliteFiles(options)

            restart_update = {"masterout":masterout,
                              "restarpos":9}

            restart_dic = writeCheckpoint(restart_update, options,
                              checkpointdata=restart_dic)
        else:
            masterout = restart_dic["masterout"]

            log = open(options.log,"a",0)
            log.write("\n %s: loaded results, you should never really be here " \
                      "unless the cleanup failed!" % timeStr())
            log.close()

    if not options.bootmaster:

        if restart_dic["restarpos"]<4:
            nboot = 1
            if options.nbootstrap is not None:
                nboot = int(options.nbootstrap)
                log = open(options.log,"a",0)
                log.write("\n %s: running in inline bootstrap mode with %i " \
                          "iterations" % (timeStr(), nboot))
                log.close()
            else:
                log = open(options.log,"a",0)
                log.write("\n %s: running in single-run mode (no bootstrap)" % timeStr())
                log.close()
            i=0

            restart_update = {"nbootstrap_i":i,
                              "nboot":nboot,
                              "restarpos":4}

            restart_dic = writeCheckpoint(restart_update, options,
                              checkpointdata=restart_dic)
        else:
            i = restart_dic["nbootstrap_i"]
            nboot = restart_dic["nboot"]

            log = open(options.log,"a",0)
            log.write("\n %s: loaded bootstrap details:" % timeStr())
            log.write("\n\t\tnumber of bootstrap runs: %i" % nboot)
            log.write("\n\t\trestarting at bootstrap run: %i" % (i+1))
            log.close()

        log = open(options.log,"a",0)
        log.write("i: %i, nboot: %i" % (i, nboot))
        log.write("restart position: %i" % restart_dic["restarpos"])
        log.close()

        DEfiles=[]
        featureset=None
        while i<nboot:

            if nboot>1:
                log = open(options.log,"a",0)
                log.write("\n\n ##############################\n" \
                          " %s: Bootstrap iteration %05i\n" \
                          " ##################################\n" \
                          "" % (timeStr(), i+1))
                log.close()

            if restart_dic["restarpos"]<5:

                log = open(options.log,"a",0)
                log.write("1")
                log.close()

                if options.nbootstrap is not None:
                    log = open(options.log,"a",0)
                    log.write("2")
                    log.close()
                    selection = rep_subselect(experiment, includelist,
                                              options, exprs_data)
                    log = open(options.log,"a",0)
                    log.write("3")
                    log.close()
                    sub_experiment = selection[0]
                    sub_all_replicates = selection[1]
                    sub_exprs_data = selection[2]

                else:
                    sub_experiment = experiment
                    sub_all_replicates = all_replicates
                    sub_exprs_data = exprs_data

                restart_update = {"sub_experiment":sub_experiment,
                                  "sub_all_replicates":sub_all_replicates,
                                  "sub_exprs_data":sub_exprs_data,
                                  "restarpos":5}

                restart_dic = writeCheckpoint(restart_update, options,
                                  checkpointdata=restart_dic)
            else:
                sub_experiment = restart_dic["sub_experiment"]
                sub_all_replicates = restart_dic["sub_all_replicates"]
                sub_exprs_data = restart_dic["sub_exprs_data"]

                log = open(options.log,"a",0)
                log.write("\n %s: loaded run sub-selection details." % timeStr())
                log.close()

            # define full set of DE features - this is needed for constructing
            # the databases later.
            if featureset is None:
                featureset = sub_exprs_data[1]

            # ==========
            # Normalize the readcount data in-script (as opposed to in algorithm)
            # if required.
            #     - allow for restart
            # ==========
            normdata = None
            if restart_dic["restarpos"]<6:
                normdata, normvals = normalize_data((sub_exprs_data[0],
                                                     sub_exprs_data[1],
                                                     sub_exprs_data[3]),
                                                    sub_experiment,
                                                    sub_all_replicates,
                                                    options
                                                    )

                restart_update = {"normdata":normdata,
                                  "normvals":normvals,
                                  "restarpos":6}

                restart_dic = writeCheckpoint(restart_update, options,
                                  checkpointdata=restart_dic)
            else:
                normdata = restart_dic["normdata"]
                normvals = restart_dic["normvals"]

                log = open(options.log,"a",0)
                log.write("\n %s: loaded normalizations" % timeStr())
                log.close()

            # ==========
            # Write the data and experiment structure to files R can use to make
            # an expressionset
            #     - allow for restart
            # ==========
            rexpfile = None
            rphefile = None
            if restart_dic["restarpos"]<7:
                rexpfile,rphefile = write_R_data_file((normdata,
                                                       sub_exprs_data[1],
                                                       sub_exprs_data[3]
                                                       ),
                                                      sub_experiment,
                                                      options
                                                      )

                restart_update = {"rexpfile":rexpfile,
                                  "rphefile":rphefile,
                                  "restarpos":7}

                restart_dic = writeCheckpoint(restart_update, options,
                                  checkpointdata=restart_dic)
            else:
                rexpfile = restart_dic["rexpfile"]
                rphefile = restart_dic["rphefile"]

                log = open(options.log,"a",0)
                log.write("\n %s: loaded R data files" % timeStr())
                log.close()

            # ==========
            # Call the appropriate R script to perform Differential Expression
            # Analysis. These scripts typically load the data and experiment
            # data and construct an expressionset object, that they then operate,
            # performing their own internal normalization if requested.
            #     - allow for restart
            # ==========
            if restart_dic["restarpos"]<8:
                
                if nboot>1:
                    fileout_details = os.path.splitext(options.outfile)
                    this_fileout = "%s_%i%s" % (fileout_details[0],
                                                i+1,
                                                fileout_details[1]
                                                )
                else:
                    this_fileout = options.outfile

                log = open(options.log,"a",0)
                log.write("\n %s: DE calls will be written to %s" \
                          "" % (timeStr(),this_fileout))
                log.close()

                # Enable running the AGNC for non-bootstrap mode if agncpath is specified
                run_AGNC=False
                if options.agncpath != None and options.nbootstrap == None:
                  run_AGNC=True

                DE_results = runRDE(rexpfile,
                                    rphefile,
                                    options,
                                    fileout=this_fileout,
                                    runAGNC=run_AGNC
                                    )

                DEfiles.append(DE_results)

                restart_update = {"bs_outfiles":DEfiles,
                                  "restarpos":8}

                restart_dic = writeCheckpoint(restart_update, options,
                                  checkpointdata=restart_dic)
            else:
                options.outfile = restart_dic["outfile"]
                DEfiles = restart_dic["bs_outfiles"]
                log = open(options.log,"a",0)
                log.write("\n %s: loaded output data" % timeStr())
                log.close()

            i+=1
            restart_update = {"restarpos":4, "nbootstrap_i":i}
            restart_dic = writeCheckpoint(restart_update, options,
                              checkpointdata=restart_dic)

    # ==========
    # Finish with housekeeping and cleanup
    # ==========

    log = open(options.log,"a",0)
    if options.bootslave:
        dbfile = results2sqlite(featureset, DEfiles, options)

        # get final fileout, munging for whether its been on the cluster or not
        try:
            masterfilename = os.path.splitext(os.path.basename(options.masterfileout))
        except AttributeError:
            options.masterfileout = options.outfile
            options.masterworkingdir = os.path.dirname(options.masterfileout)
            masterfilename = os.path.splitext(os.path.basename(options.masterfileout))
                
        masteroutfile = "%s/%s.db" % (options.masterworkingdir,
                                      masterfilename[0])
                
        log.write("\n %s: Copying database file (%s), to master working " \
                  "dir as file (%s)...\n" % (timeStr(),
                                     dbfile,
                                     masteroutfile
                                     )
                  )

        if os.path.exists(masteroutfile):
            log.write("\tremoving pre-existing version of %s... \n" \
                      "" % masteroutfile
                      )
            os.remove(masteroutfile)

        if os.path.exists(dbfile):
            shutil.copy(dbfile, masteroutfile)
            ex = os.path.exists(masteroutfile)
            if not ex:
                log.write(" %s: ERROR: Something went wrong copying the " \
                          "database file to the master directory. Get help!")
        else:
            log.write("help! I've lost the %s file...\n" % dbfile)
            log.write("\tcwd: %s\n" % os.getcwd())
            for dfilename in os.getcwd():
                log.write("\t\t%s\n" % dfilename)
            log.write("\ttmpdir: %s\n" % options.tmpdir)
            for dfilename in options.tmpdir:
                log.write("\t\t%s\n" % dfilename)
    
    elif options.nbootstrap is None:
        log.write("\n %s: Moving DE call output (%s) to final output (%s)" \
                  % (timeStr(), DEfiles[0], options.outfile))
        shutil.move(DEfiles[0], options.outfile)
    
    if options.keep_tmpdir:
        log.write("\n %s: --keep_tmpdir option set. Temporary files can be " \
                  "found here: %s" % (timeStr(),options.tmpdir)
                  )
        log.write("\n All finished. Enjoy... :)\n")
    else:
        createNewTempdir(options, remove=True, create=False)
        log.write("\n %s: Finished. Enjoy... :)\n" % timeStr())

    log.close()

##
## TODO - in concatSqlliteFiles() need a proper check that the correct number of
##        db files were found. No checks are done at all ATM.
##      - Other normalisations need to be coded up.
##      - Calculate R-script runtime and store.
##      - Restart does not work for once passed the the point the bootstrap is
##        called if options.nbootstrap is None as the "nbootstrap_i" is the
##        same as the "nboot" value.
##
## FIXED
##       - fixed bug with final output file not being generated if run in 
##         non-bootstrap mode
##       - updates species default. fixed bug with not excluding files correctly
##         if using precounts
##       - updated the way tempfiles and dirs work so that when -v is specified
##         all the tmpfiles are in tmpdirs in the local tmpdir rather than on
##         individual cluster nodes. Makes it much easier to debug shit.
##       - updated documentation and usage strings for all functions and the
##         script as a whole, including a description of the output structure.
##       - low priority: more complex job status monitoring for cluster jobs.
##       - seems like we can't depend on SGE exit codes for 'SUCCESS' as some
##           jobs core-dump and exit succesfully. Add checking for expected files.
##       - add_gene_name_column.pl is unnecessarily run for every
##           group_by_gene.pl output. Turn this off.
##       - tmpdir option now automatically resolves relative paths.
##       - the drmaa session now gets passed the current working directory and
##           the -l local_disk=10G option to make sure there is free space on
##           the node before it runs.
##       - The default for clustq is now hard-coded to be 64bit-pri (this will
##           need to be changed is running on a different cluster) and is
##           submitted as a parameter for the drmaa runs every time.
##       - Added PERLLIB, PERL5LIB, R_HOME, PATH, DRMAA_LIBRARY_PATH, and
##           LD_LIBRARY_PATH to the environment variables logged by the code
##       - Added R version info to log header
##       - Pre-existing tmpdirs are not removed when a new run starts, if found.
##       - Updated group_by_gene.pl to call htseq-count in quiet mode. Hoped for
##           speed up of runtime. No real benefit. But her, no more 2.5GB of
##           logs...
##       - uses localtime() instead of gmtime()
##       - capture version info from R script
##       - added explicit perlpath specification option and version logging
##           perl scripts are now run using this specific perl version
##       - fixed a major bug that deleted root bam files --- oops!
##       - reformated docstrings.
##       - Added sanity checking to make sure the bootstrap selects unique
##         replicates between conditions. This should be irrellavant unless the
##         false discovery rate to tools is being tested by having two identical
##         conditions (13/09/17)
##
