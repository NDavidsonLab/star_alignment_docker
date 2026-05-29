#!/bin/bash

script_dir=$(dirname "$(realpath "$0")")
cd "${script_dir}/../"

# ---

# inputs
genome_dir="./genome/HG38/STAR"
data_dir="./data"

# output
outdir="./alignments"
mkdir -p ${outdir}
mkdir -p ${outdir}/fastqc
mkdir -p ${outdir}/trimmed
mkdir -p ${outdir}/bams

# Runs FastQC quality control on a pair of FASTQ files (or a BAM file) for a given sample.
# Arguments:
#   $1 - prefix: sample name/identifier, used to name the output subdirectory
#   $2 - fastq1: path to the first FASTQ file (or BAM file)
#   $3 - fastq2: path to the second FASTQ file (optional; if absent, input is assumed to be a BAM)
run_fastqc () {
    local prefix=$1
    local fastq1=$2
    local fastq2=$3

    # Log the inputs and output destination for this run
    echo "* in run_fastqc(), for prefix ${prefix}:"
    echo "- Processing fastq1: ${fastq1}"
    echo "- Processing fastq2: ${fastq2}"
    echo "- Output to ${outdir}/fastqc/${prefix}"
    echo ""

    # fastq1 is required; exit early if not provided
    if [ -z "${fastq1}" ]; then
        echo "missing fastq1"
        return
    fi

    # fastq2 is optional; a missing second file is treated as BAM input
    if [ -z "${fastq2}" ]; then
        echo "missing fastq2, assuming BAM"
    fi

    # Create the per-sample FastQC output directory if it doesn't already exist
    if [ ! -d ${outdir}/fastqc/${prefix} ]; then
        mkdir ${outdir}/fastqc/${prefix}
    fi

    # Run FastQC with 12 threads, writing results to the sample-specific output directory
    echo "RUN"
    fastqc -t 12 -o ${outdir}/fastqc/${prefix} ${fastq1} ${fastq2}
}

# Runs fastp adapter trimming and quality filtering on a pair of paired-end FASTQ files.
# Both fastq1 and fastq2 are required; the function exits early if either is missing.
# Trimmed output files are written to ${outdir}/trimmed/ with _R1.trimmed and _R2.trimmed suffixes.
# Arguments:
#   $1 - prefix: sample name/identifier, used to name the output files
#   $2 - fastq1: path to the R1 (forward) FASTQ file
#   $3 - fastq2: path to the R2 (reverse) FASTQ file
run_trimm () {
    local prefix=$1
    local fastq1=$2
    local fastq2=$3

    # Log the inputs and output destination for this run
    echo "* in run_trimm(), for prefix ${prefix}:"
    echo "- Processing fastq1: ${fastq1}"
    echo "- Processing fastq2: ${fastq2}"
    echo "- Output to ${outdir}/trimmed/${prefix}"
    echo ""

    # Both FASTQ files are required for paired-end trimming; exit early if either is missing
    if [ -z "${fastq1}" ]; then
        echo "missing fastq1"
        return
    fi

    if [ -z "${fastq2}" ]; then
        echo "missing fastq2"
        return
    fi

    # Create the trimmed output directory if it doesn't already exist
    if [ ! -d ${outdir}/trimmed/ ]; then
        mkdir ${outdir}/trimmed
    fi

    # Run fastp with paired-end adapter auto-detection, poly-G tail trimming,
    # sliding-window quality trimming from the tail end, and 8 threads
    echo "RUN"
    fastp \
        -i ${fastq1} -I ${fastq2} \
        -o ${outdir}/trimmed/${prefix}_R1.trimmed.fastq.gz -O ${outdir}/trimmed/${prefix}_R2.trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --cut_tail \
        --thread 8
}


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

    if [ -f ${outdir}/bams/${prefix}/${prefix}Aligned.sortedByCoord.out.bam ]; then
        echo "${outdir}/bams/${prefix}/${prefix}Aligned.sortedByCoord.out.bam already exists, skipping"
        return
    fi

    if [ ! -d ${outdir}/bams/${prefix}/ ]; then
        mkdir ${outdir}/bams/${prefix}
    fi

    echo "RUN"

    STAR \
        --runThreadN 12 \
        --genomeDir ${genome_dir} \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${fastq1} ${fastq2} \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI NM MD AS XS \
        --outSAMunmapped Within \
        --outSAMheaderHD @HD VN:1.4 \
        --outFilterMultimapNmax 20 \
        --outFilterMultimapScoreRange 1 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --alignSJDBoverhangMin 1 \
        --sjdbOverhang 100 \
        --sjdbScore 2 \
        --outFileNamePrefix ${outdir}/bams/${prefix}/${prefix} \
        --readFilesCommand zcat \
        --quantMode GeneCounts TranscriptomeSAM

    echo ""
    echo "* Sorting BAM with samtools..."
    samtools sort \
        -m 4G \
        -@ 12 \
        -o ${outdir}/bams/${prefix}/${prefix}Aligned.sortedByCoord.out.bam \
        ${outdir}/bams/${prefix}/${prefix}Aligned.out.bam

    echo "* Indexing sorted BAM..."
    samtools index ${outdir}/bams/${prefix}/${prefix}Aligned.sortedByCoord.out.bam

    echo "* Removing unsorted BAM..."
    rm ${outdir}/bams/${prefix}/${prefix}Aligned.out.bam

    echo ""
}

##### RUN

# --- Step 1: Trim adapters and low-quality bases ---
# Collect all FASTQ files in the data directory and extract sample prefixes
# by cutting on '_' and keeping fields 1-3 (e.g. /path/to/SAMPLE_S1_L001)
prefixes=( $( ls ${data_dir}/*/*.fq.gz | cut -d'_' -f1 ) )
# Deduplicate the prefix list (each sample appears once per R1/R2 file otherwise)
prefixes_uniq=($(echo "${prefixes[@]}" | tr ' ' '\n' | sort -u))

for prefix in ${prefixes_uniq[@]}; do
    echo "QCing ${prefix} paired files: "
    echo ${prefix}*
    echo ""

    # Optional first-pass FastQC on raw reads (currently commented out)
    echo "First pass fastqc"
    run_fastqc $( basename ${prefix} ) ${prefix}*

    # Trim adapters and low-quality bases with fastp
    echo ""
    echo "Trimming"
    run_trimm $( basename ${prefix} ) ${prefix}*

    echo ""
done


# --- Step 2: Align trimmed reads to the reference genome ---
# Collect trimmed FASTQ files and extract sample prefixes the same way
prefixes=( $( ls ${outdir}/trimmed/*.fastq.gz | cut -d'_' -f1 ) )
# Deduplicate so each sample is aligned once
prefixes_uniq=($(echo "${prefixes[@]}" | tr ' ' '\n' | sort -u))

for prefix in ${prefixes_uniq[@]}; do
    echo "Aligning ${prefix} paired files: "
    echo ${prefix}*
    echo ""

    # Align with STAR; passes both trimmed R1 and R2 via glob expansion
    echo "Aligning"
    run_star $( basename ${prefix} ) ${prefix}*

    echo ""
done


# --- Step 3: Run final FastQC on aligned BAM files ---
# Collect all BAM files across per-sample subdirectories and extract the path
# prefix up to the 'A' in 'Aligned' (so the glob ${prefix}Aligned.sortedByCoord*
# resolves to the correct output BAM)
prefixes=( $( ls ${outdir}/bams/*/*.bam | cut -d'A' -f1 ) )
prefixes_uniq=($(echo "${prefixes[@]}" | tr ' ' '\n' | sort -u))

for prefix in ${prefixes_uniq[@]}; do
    echo "Final QCing ${prefix} files: "
    echo ${prefix}Aligned.sortedByCoord*
    echo ""

    # Run FastQC on the coordinate-sorted BAM to assess post-alignment quality
    echo "Final pass fastqc"
    echo $( basename ${prefix} )
    run_fastqc $( basename ${prefix} ) ${prefix}Aligned.sortedByCoord*

    echo ""
done


# --- Step 4: Aggregate QC reports with MultiQC ---
# Collect all FastQC and other QC outputs from the output directory into a
# single interactive report, written to ${outdir}/multiqc_report
mkdir -p ${outdir}/multiqc_report
multiqc ${outdir} -o ${outdir}/multiqc_report