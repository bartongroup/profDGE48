#!/bin/bash
# set you email for grid engine job failure notifications
myemail="p.schofield@dundee.ac.uk"
# you store your grid engine logs
logs="${LOGS}"
# your copy of the rna-seq code is
rnaseq="${HOME}/Projects/rna_seq/"
# your bam files are
datadir="/cluster/gjb_lab/cdr/GRNAseq/analysis/\
false_positive_rate_testing/conditions/"
# annotation file
annofile="${rnaseq}Annotations/Saccharomyces_cerevisiae.EF4.68.gtf"
# temporary directory
tmpdir="${rnaseq}NOBACK/nospikes_fdr/"
# code directory
codedir="${rnaseq}DE_Tools/"
# exlude list
exlist="${rnaseq}nonsvn/exclude_nospikes.lst"
# samtools path
samtoolspath="/sw/opt/samtools-0.1.18/samtools"

# method comes from command line
meth=${1}
START=$(eval echo "${2}")
END=$(eval echo "${3}")
norm=${4}
normid=$(echo ${norm} | cut --output-delimiter="_" -d ":" -f1,2 )  

# Run bootstraps for default normalization
for i in $(eval echo "{${START}..${END}}")
  do 
		qsub -l ram=1G -l urna=1 -o ${logs} -e ${logs} -V -b y \
      python ${codedir}generic_wrapper.py \
				-d ${datadir} \
				-a ${annofile} \
        -o ${tmpdir}${meth}_${i}reps_${normid}.db  \
        -l ${tmpdir}${meth}_${i}reps_${normid}.log \
        -r ${codedir}${meth}.R \
				-k ${i} -b 100 \
				-e ${exlist} \
        --tmpdir ${tmpdir}${meth}_${i}reps_${normid} \
        --samtoolspath ${samtoolspath} \
        --gbgfile ${codedir}group_by_gene.pl \
        --agncfile ${codedir}add_gene_name_column.pl \
				--keep-tmpdir --bootstrapmaster --verbose \
				--norm ${norm} --precounts 
  done

