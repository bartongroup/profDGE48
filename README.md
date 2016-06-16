# Profiling Differential Gene Expression with a 48 replicate experiment

This is the Github repository for the collection of scripts developed by 
members of [The Barton Group] (http://www.compbio.dundee.ac.uk) at [The 
University of Dundee] (http://www.dundee.ac.uk) to process and analyse the 
data from a 48-replicate RNA-seq experiment conducted specifically to test the
underlying assumptions and performance of popular RNA-seq Differential Gene 
Expression (DGE) tools. A full description of the experiment is given in 
[Schurch et. al. 2015] (https://doi.org/10.1261/rna.053959.115).

## Installation

The scripts in this collection are written in a variety of languages including 
[perl](https://www.perl.org/), [python](https://docs.python.org/2/) and 
[R](https://www.r-project.org/). None of the scripts in the package require 
any specific installation however there are a considerable number of dependencies
for the full suite of scripts including (in some cases) access to a DRMAA-enabled 
cluster and they have only been tested with specific versions of the languages and 
dependencies. For more information on the specific requirements see the main 
documentation in the codebase (*Doc/html*).

## Getting started

This is a very brief overview of how to get started using the code. For a detailed 
walk-through please refer to the *Getting Started* section of the codebase documentation
(*Doc/html/getting_started.html*). Before running any of the scripts you will need to 
clone the repository, set up the appropriate environmental variables *PERL5LIB*, 
*PYTHONPATH*, and *PATH*, and you will need some data. 

### Data

The codebase includes pre-processed intermediate level data or the original
48-replicate experiment (in the *Preprocessed_data* folder). Alternatively the raw 
fastq data for the experiment can be obtained from the 
[European Nucleotide Archive](http://www.ebi.ac.uk/ena/data/view/PRJEB5348)
or any sufficiently replicated fastq data can be used. If you are using fastq data 
the data will need to be aligned to the relevant genome and provided as indexed 
*.bam* files (the genome sequence used for the 48-replicate experiment is available 
in the Annotations folder of the codebase).

### Gold-standard Differential Gene Expression (DGE) results

The performance of each DGE tool as a function of replicate number and expression 
fold-change is evaluated by comparing the DGE results from sub-selections of the 
replicate data against a 'gold standard' set of DGE results calculated for each tool 
using the full set of (clean) replicates. The codebase includes pre-computed gold
standards for each of the tools tested in the 48 replicate experiment (in the 
*Preprocessed_data/full_rep_gold_standards* folder). Alternatively, gold standards 
can be computed using the *Bootstrapping/generic_wrapper.py* script with details of 
which tool to use. A simple example command-line for generating a gold standard for 
*edgeR* might be something line:

> Bootstrapping/generic_wrapper.py -r DE_tool_scripts/edgeR.R -d data_dir -a annotation.gtt -o golds/edgeR_gold.txt -l edgeR_gold.log

### Bootstrap DGE results

The DGE from bootstrapped sub-selections of replicates is computed in a similar 
fashion to the gold standards, however we now additionally specify a sqlite output format,
a number of replicates to select and a number of bootstraps to perform. For the analysis 
in [Schurch et. al. 2015] (http://arxiv.org/abs/1505.02017) use 100 bootstraps and replicate
sub-selection from 2..40. Again, an example for *edgeR* with 10 bootstraps os sub-selections 
with 3 replicates in each condition might be something like:

> Bootstrapping/generic_wrapper.py -r DE_tool_scripts/edgeR.R -d data_dir -a annotation.gtf -b 10 -k 3 -o edgeR/edgeR_k03.db -l edgeR_k03.log

The output sqlite database contains the results from and log information for all the
individual bootstrap calls. These can then be compared directly with the gold standard
results. A variety of statistics and plotting scripts are available the codebase for this. 
Please see the *Getting Started* section of the codebase documentation
(*Doc/html/getting_started.html*) for a simple example of how to plot a tools performance
or see the individual documentation for the individual scripts 
(*Docs/html/package_inventory.html*) for more details.

## Contact information

For further information or assistance with this repository please contact one of:

* Nick Schurch: <nschurch@dundee.ac.uk>
* Marek Gierlinski: <mgierlinski@dundee.ac.uk>
* Pieta Schofield: <pschofield@dundee.ac.uk>
* Chris Cole: <ccole@dundee.ac.uk>
* Geoff Barton: <gjbarton@dundee.ac.uk>

