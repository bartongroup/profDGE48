#!/bin/bash
# set you email for grid engine job failure notifications
myemail="p.schofield@dundee.ac.uk"
# you store your grid engine logs
logs=${LOGS}
# your copy of the rna-seq code is
rnaseq="/cluster/gjb_lab/pschofield/Projects/grna_dag/Github_cleanup/"
# your bam files are
datadir="/cluster/gjb_lab/cdr/GRNAseq/mapping/genome_biol_reps"
# annotation file
annofile="${rnaseq}Annotations/Saccharomyces_cerevisiae.EF4.68.gtf"
# temporary directory
tmpdir="${HOME}/scratch/NOBACK/"
# code directory
codedir="${rnaseq}Bootstrapping/"
# tool code directory
tooldir="${rnaseq}DE_tool_scripts/"
# exlude list
exlist=""
# samtools path
samtoolspath="/sw/opt/samtools-0.1.18/samtools"
# rpath
rpath="/sw/opt/R/3.2.2/bin/Rscript"

# method comes from command line
meth=${1}
START=$(eval echo "${2}")
END=$(eval echo "${3}")
norm=${4}
normid=$(echo ${norm} | cut --output-delimiter="_" -d ":" -f1,2 )

# Run bootstraps for default normalization
module load python/2.7.3
module load perl

export PERL5LIB="/sw/lib/perl5.10.1/lib/perl5:${PERL5LIB}"
export PYTHONPATH="${rnaseq}:${PYTHONPATH}"

# Run bootstraps for default normalization
for i in $(eval echo "{${START}..${END}}")
  do
    methdir="${tmpdir}/${meth}/"
    rtmpdir="${tmpdir}/${meth}/${meth}_${i}reps_${normid}/"
    mkdir -p ${rtmpdir}
    qsub -l ram=1G -l urna=1 -o ${logs} -e ${logs} -V -q c6100.q -b y \
      python ${codedir}generic_wrapper.py \
				-d ${datadir} \
				-a ${annofile} \
        -o ${methdir}${meth}_${i}reps_${normid}.db  \
        -l ${methdir}${meth}_${i}reps_${normid}.log \
        -r ${tooldir}${meth}.R \
				-k ${i} -b 100 \
				--Rpath ${rpath} \
        --tmpdir ${rtmpdir} \
        --samtoolspath ${samtoolspath} \
        --gbgfile ${codedir}group_by_gene.pl \
        --agncfile ${codedir}add_gene_name_column.pl \
				--keep-tmpdir --bootstrapmaster --verbose \
				--norm ${norm} \
				--species s_cerevisiae \
				--precounts --savecounts
  done
