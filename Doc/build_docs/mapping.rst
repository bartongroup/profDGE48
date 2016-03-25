.. mapping:

.. _mapping:

*******************************
Mapping reads and summary stats
*******************************

.. build_pileup_perldoc:

.. _build_pileup_perldoc:

***************
build_pileup.pl
***************

Build aggregated pileup from all BAM files::

	build_pileup.pl -datapath=./mydata

:param: `-datapath|d` - Path to the data directory of the experiment. Each sub-directory in this will be treated as a separate condition in the experiment (with the name of the condition matching the name of the directory), and each .bam file in each directory is a replicate for that condition.

:param: `-outpath|p` - Path to the output directory. 

:param: `-conditions|c` - Comma-separated list of conditions. (Default: *WT,Snf2*)

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. concat_fastq_perldoc:

.. _concat_fastq_perldoc:

***************
concat_fastq.pl
***************

This script takes the multilane, technical replicate data (as produced by 
GenePool) and combines it into biological replicate fastq files. It assumes 
the source directory has 96 subdirs (numbered 1-96) and that the first 48 
are 'Snf2' samples and the second 48 are 'WT' samples. It also assumes that 
in each of the number subdirs there are seven technical replicate files, 
which need combining into one, biological, replicate.

Output files are put in the destination directory under 'Snf2' and 'WT' 
subdirs. Filenames are derived from the biological replicate number, 
sample and MID. e.g. *WT_rep06_MID02_allLanes.fastq*::

	concat_fastq.pl --src-path <path> --dest-path <path> [--version] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

:param: `--src-path` - Source directory with all replicate data. Must have numeric subdirs with fastq files.

:param: `--dest-path` - Destination directory for saving concatenated fastq files

:param: `--version` - Report version info and exit

:param: `--verbose|--no-verbose` - Toggle verbosity. [default:none]

:param: `--debug|--no-debug` - Toggle debugging output. [default:none]

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.
				
.. count_fastq_reads_lanes_perldoc:

.. _count_fastq_reads_lanes_perldoc:

**************************
count_fastq_reads_lanes.pl
**************************

Counts reads in each condition, replicate and lane. Creates a file 
*all_readcount_lane.dat* with counts. Requires directory with subdirs, one 
per biological sample, containing lean FASTQ files. There are some numbers and 
names hardcoded for our data, so this simple script needs redoing for other 
data sets::

	count_fastq_reads_lanes.pl -dir=/dir/with/fastq/files
      
:param: `-dir` - Directory with FASTQ files. File names should be in a specific format, so this very simple script really works only with our experiment.

.. count_fastq_reads_perldoc:

.. _count_fastq_reads_perldoc:

********************
count_fastq_reads.pl
********************

Counts reads in all fastq files and store results in output files 
*<cond>_readcount.dat*.There are some numbers and names hardcoded for our 
data, so this simple script needs redoing for other data sets::

	count_fastq_reads.pl -dir=/dir/with/fastq/files
      
:param: `-dir` - Directory with FASTQ files. File names should be in a specific format, so this very simple script really works only with our experiment.

.. make_bamlinks_perldoc:

.. _make_bamlinks_perldoc:

****************
make_bamlinks.pl
****************

Creates symbolic links to .bam files in subdirectories. 
:ref:`run_alignements.pl <run_alignments_perldoc>` and 
:ref:`run_biological_alignements.pl <run_biological_alignments_perldoc>` 
scripts create individual directories for each replicate under one top 
(condition) directory. Each directory contain C<accepted_hits.bam> file. 
However, further analysis requires that all bam files are in one directory. 
This script creates symbolic links from the top directory to individual bam 
files.

  make_bamlinks.pl -topdir=/cluster/gjb_lab/cdr/GRNAseq/mapping/genome/WT

:param: `-topdir` - Top directory containing bam subdirectories. Typically, this contains one condition.

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. pileup_var_perldoc:
		
.. _pileup_var_perldoc:

*************
pileup_var.pl
*************

Create files with pileup statistic for each gene: one file with mean and one 
with standard deviation. There is one gene per line, and each line contains 
gene name and 48 numbers for all replicates. Requires pileup files for each 
chromosome, created with :ref:`build_pileup.pl <build_pileup_perldoc>`.

  pileup_var.pl

:param: `-pileupdir` - The directory with pileups (created with :ref:`build_pileup.pl <build_pileup_perldoc>`). If not defined, the default value from *defs.dat* file will be used.

:param: `-outdir` - The directory to save results. If not specified, it will the pileup directory.

:param: `-gff` - The GTF file to be used. The default value is *Saccharomyces_cerevisiae.EF4.68.gtf*. 

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. run_alignments_perldoc:

.. _run_alignments_perldoc:

*****************
run_alignments.pl
*****************

This scripts takes all 96 replicates from the yeast snf2 mutant experiment 
and aligns the raw Fastq reads to the genome with `Tophat <https://ccb.jhu.edu/software/tophat/index.shtml>`_. 
The *--src*  directory must have 96 subdirs numbered 1-96, where 1-48 are the snf2 
mutants and 49-96 are the wild-types. 

The *--dest* directory must have the *Snf2/* and *WT/* subdirs where the 
`Tophat <https://ccb.jhu.edu/software/tophat/index.shtml>`_ runs will be saved::

	run_alignments.pl --src <dir> --dest <dir> [--genome <pathprefix>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

:param: `--src` - Source path of subdirs with fastq files 

:param: `--dest` - Destination path for output.

:param: `--genome` - Path and to bowtie2 genome index [default: /db/bowtie2/Scerevisiae68_ERCC92]

:param: `--verbose|--no-verbose` - Toggle verbosity. [default:none]

:param: `--debug|--no-debug` - Toggle debugging output. [default:none]

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. run_biological_alignments_perldoc:

.. _run_biological_alignments_perldoc:

****************************
run_biological_alignments.pl
****************************

Given a set of combined biological replicate fastq files run 
`Tophat <https://ccb.jhu.edu/software/tophat/index.shtml>`_ on them. Use this 
to run `Tophat <https://ccb.jhu.edu/software/tophat/index.shtml>`_ over the set 
of fastq files where technical replicates have been combined into biological 
replicates. For this experiment it means 48 replicates each of WT and 
Snfs2 samples::

	run_biological_alignments.pl --src-path <path> --dest-path <path> [--verbose|--no-verbose] [--version] [--debug|--no-debug] [--man] [--help]

:param: `--src-path` - Path containing 'WT' and 'Snf2' subdirs with fastq files.

:param: `--dest-path` - Path containing 'WT' and 'Snf2' subdirs for saving tophat output.

:param: `--version` - Report version information and exit

:param: `--verbose|--no-verbose` - Toggle verbosity. [default:none]

:param: `--debug|--no-debug` - Toggle debugging output. [default:none]

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. sort_and_index_bams_autodocs:

.. _sort_and_index_bams_autodocs:

**********************
sort_and_index_bams.py
**********************

.. automodule:: sort_and_index_bams

=========
Functions
=========

.. autofunction:: sortOrIndexBam
.. autofunction:: sortOrIndexBams
		
.. summary_stats_perldoc:

.. _summary_stats_perldoc:

****************
summary_stats.pl
****************

Generate summary statistics of alignment infomation. This script generates 
simple statistics on the outcome of Tophat alignments. It expects that the 
given *--path* includes subdirectories generated by tophat which include 
*accepted_hits.bam* and *unmapped.bam* files in them. Given that assumption 
the script with then run '`samtools <http://www.htslib.org/>`_ *flagstat*' on 
the BAM files to determine the number of mapped, unmapped and failed reads for 
each run. These are then reported in the output file::

	summary_stats.pl --path <path> [--samtools <path>] [--threads <num>] [--out <file>] [--version] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

:param: `--path` - Path to tophat output dirs.

:param: `--samtools` - Path to samtools executeable [default: /sw/samtools-0.1.18/samtools]

:param: `--threads` - Number of threads to use [default: 1]

:param: `--out` - Output filename. [default: out.tsv]

:param: `--version` - Output version info and exit.

:param: `--verbose|--no-verbose` - Toggle verbosity. [default: verbose]

:param: `--debug|--no-debug` - Toggle debugging output. [default: no-debug]

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

.. unique_fastq_reads_perldoc:

.. _unique_fastq_reads_perldoc:

*********************
unique_fastq_reads.pl
*********************

Takes a Fastq file and writes out only unique examples of reads therein. 
Keeps header and quality information for the first instance of each read::

	unique_fastq_reads.pl --in <file> [--out <file>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

:param: `--in` - Input fastq file (can be gzipped).

:param: `--out` - Output filename. [default: unique.fastq]

:param: `--verbose|--no-verbose` - Toggle verbosity. [default:none]

:param: `--debug|--no-debug` - Toggle debugging output. [default:none]

:param: `--help` - Brief help.

:param: `--man` - Full manpage of program.

..	toctree::
	:maxdepth: 1
