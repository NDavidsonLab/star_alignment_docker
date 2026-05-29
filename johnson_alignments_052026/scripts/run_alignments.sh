#!/bin/bash

script_dir=$(dirname "$(realpath "$0")")
cd "${script_dir}/.."

# ---

# inputs
genome_dir="./genome/M37/STAR"
data_dir="./data"

# output
outdir="./alignments"
mkdir -p ${outdir}

run_star () {
    local prefix=$1
    local fastq1=$2
    local fastq2=$3

    echo "* in run_star(), for prefix ${prefix}:"
    echo "- Processing fastq1: ${fastq1}"
    echo "- Processing fastq2: ${fastq2}"
    echo "- Output to ${outdir}/${prefix}"
    echo ""

    if [ -z "${fastq1}" ]; then
        echo "missing fastq1"
        return
    fi

    if [ -z "${fastq2}" ]; then
        echo "missing fastq2"
        return
    fi

    if [ -f ${outdir}/${prefix}Aligned.sortedByCoord.out.bam ]; then
        echo "${outdir}/${prefix}Aligned.sortedByCoord.out.bam already exists, skipping"
        return
    fi

    echo "RUN"

    STAR \
        --runThreadN 12 \
        --genomeDir ${genome_dir} \
        --genomeLoad NoSharedMemory  \
        --readFilesIn ${fastq1} ${fastq2} \
        --outSAMtype BAM  SortedByCoordinate \
        --outSAMstrandField intronMotif  \
        --outSAMattributes NH  HI  NM  MD  AS  XS \
        --outSAMunmapped Within \
        --outSAMheaderHD @HD  VN:1.4 \
        --outFilterMultimapNmax 20  \
        --outFilterMultimapScoreRange 1 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --alignSJDBoverhangMin 1 \
        --sjdbOverhang 100 \
        --sjdbScore 2 \
        --outFileNamePrefix ${outdir}/${prefix} \
        --limitBAMsortRAM 30000000000 \
        --readFilesCommand zcat

    echo ""
}

prefixes=( $( ls ${data_dir}/*.fastq.gz | cut -d'_' -f1-4 ) )

for prefix in ${prefixes[@]}; do
    echo "Processing ${prefix} paired files: "
    echo ${prefix}*
    echo ""
    run_star $( basename ${prefix} ) ${prefix}*
    echo ""
done
