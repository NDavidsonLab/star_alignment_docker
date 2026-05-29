#!/bin/bash

script_dir=$(dirname "$(realpath "$0")")
cd "${script_dir}/../"

# ---

# inputs
genome_dir="./genome/M37/STAR"
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


# Aligns paired-end FASTQ files to a reference genome using STAR and outputs a
# coordinate-sorted BAM file along with gene counts and transcriptome-aligned reads.
# Skips alignment if the output BAM already exists.
# Arguments:
#   $1 - prefix: sample name/identifier, used to name the output directory and files
#   $2 - fastq1: path to the R1 (forward) FASTQ file (gzipped)
#   $3 - fastq2: path to the R2 (reverse) FASTQ file (gzipped)
run_star () {
    local prefix=$1
    local fastq1=$2
    local fastq2=$3

    # Log the inputs and output destination for this run
    echo "* in run_star(), for prefix ${prefix}:"
    echo "- Processing fastq1: ${fastq1}"
    echo "- Processing fastq2: ${fastq2}"
    echo "- Output to ${outdir}/${prefix}"
    echo ""

    # Both FASTQ files are required for paired-end alignment; exit early if either is missing
    if [ -z "${fastq1}" ]; then
        echo "missing fastq1"
        return
    fi

    if [ -z "${fastq2}" ]; then
        echo "missing fastq2"
        return
    fi

    # Skip alignment if the sorted BAM output already exists for this sample
    if [ -f ${outdir}/bams/${prefix}${prefix}/Aligned.sortedByCoord.out.bam ]; then
        echo "${outdir}/bams/${prefix}/${prefix}Aligned.sortedByCoord.out.bam already exists, skipping"
        return
    fi

    # Create the per-sample BAM output directory if it doesn't already exist
    if [ ! -d ${outdir}/bams/${prefix}/ ]; then
        mkdir ${outdir}/bams/${prefix}
    fi

    echo "RUN"

    # --runThreadN 12                    : use 12 threads
    # --genomeDir                        : path to STAR genome index
    # --genomeLoad NoSharedMemory        : load genome into memory fresh each run
    # --readFilesIn                      : paired-end input FASTQs
    # --outSAMtype BAM SortedByCoordinate: output coordinate-sorted BAM
    # --outSAMstrandField intronMotif    : add strand field based on intron motif (for unstranded libraries)
    # --outSAMattributes                 : include standard alignment attributes in BAM
    # --outSAMunmapped Within            : include unmapped reads in the BAM output
    # --outSAMheaderHD                   : set SAM header version
    # --outFilterMultimapNmax 20         : discard reads mapping to more than 20 locations
    # --outFilterMultimapScoreRange 1    : allow multimappers within score range of 1 of the best
    # --outFilterScoreMinOverLread 0.33  : minimum alignment score as fraction of read length
    # --outFilterMatchNminOverLread 0.33 : minimum matched bases as fraction of read length
    # --outFilterMismatchNmax 10         : maximum number of mismatches per read
    # --alignIntronMax 500000            : maximum intron size
    # --alignMatesGapMax 1000000         : maximum gap between paired-end mates
    # --alignSJDBoverhangMin 1           : minimum overhang for annotated splice junctions
    # --sjdbOverhang 100                 : overhang length for splice junction database (typically read length - 1)
    # --sjdbScore 2                      : extra alignment score for reads spanning annotated junctions
    # --outFileNamePrefix                : prefix for all output files
    # --limitBAMsortRAM 30000000000      : RAM limit (30 GB) for BAM sorting
    # --readFilesCommand zcat            : decompress gzipped FASTQs on the fly
    # --quantMode GeneCounts TranscriptomeSAM : output gene counts table and transcriptome-aligned BAM
    STAR \
        --runThreadN 12 \
        --genomeDir ${genome_dir} \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${fastq1} ${fastq2} \
        --outSAMtype BAM SortedByCoordinate \
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
        --limitBAMsortRAM 30000000000 \
        --readFilesCommand zcat \
        --quantMode GeneCounts TranscriptomeSAM
    echo ""
}

##### RUN

# --- Step 1: Trim adapters and low-quality bases ---
# Collect all FASTQ files in the data directory and extract sample prefixes
# by cutting on '_' and keeping fields 1-3 (e.g. /path/to/SAMPLE_S1_L001)
prefixes=( $( ls ${data_dir}/*.fastq.gz | cut -d'_' -f1-3 ) )
# Deduplicate the prefix list (each sample appears once per R1/R2 file otherwise)
prefixes_uniq=($(echo "${prefixes[@]}" | tr ' ' '\n' | sort -u))

for prefix in ${prefixes_uniq[@]}; do
    echo "QCing ${prefix} paired files: "
    echo ${prefix}*
    echo ""

    # Optional first-pass FastQC on raw reads (currently commented out)
    echo "First pass fastqc"
    #run_fastqc $( basename ${prefix} ) ${prefix}*

    # Trim adapters and low-quality bases with fastp
    echo ""
    echo "Trimming"
    run_trimm $( basename ${prefix} ) ${prefix}*

    echo ""
done


# --- Step 2: Align trimmed reads to the reference genome ---
# Collect trimmed FASTQ files and extract sample prefixes the same way
prefixes=( $( ls ${outdir}/trimmed/*.fastq.gz | cut -d'_' -f1-3 ) )
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