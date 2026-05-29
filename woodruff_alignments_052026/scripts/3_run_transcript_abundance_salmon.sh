#!/bin/bash

# Runs Salmon in alignment-based mode on STAR transcriptome-aligned BAMs to
# generate transcript-level abundance estimates (TPM and expected counts).
# Requires:
#   - STAR run with --quantMode TranscriptomeSAM to produce Aligned.toTranscriptome.out.bam
#   - A transcriptome FASTA file for the Salmon index (or a pre-built index)
#   - Salmon installed and on PATH
#
# Output: ${outdir}/salmon/${sample_name}/
#   - quant.sf         transcript-level TPM and counts (import into tximport)
#   - quant.genes.sf   gene-level summary
#   - aux_info/        auxiliary files including model parameters and QC stats

# --- configuration ---
script_dir=$(dirname "$(realpath "$0")")
cd "${script_dir}/../"

indir="./alignments"

outdir="./abundances"
mkdir -p ${outdir}

genome_dir="./genome/M37"
TRANSCRIPTOME_FASTA="${genome_dir}/generated_transcriptome_for_salmon.fa"   # FASTA of transcript sequences
THREADS=12
SALMON_OUTDIR="${outdir}/salmon"
mkdir -p ${outdir}/salmon


# --- generate transcriptome FASTA from genome and GTF ---
# This ensures transcript IDs match those in the STAR transcriptome BAM,
# since both are derived from the same genome FASTA and GTF

GTF_FILE="${genome_dir}/gencode.vM37.primary_assembly.annotation.gtf"
GENOME_FASTA="${genome_dir}/GRCm39.primary_assembly.genome.fa"

if [ ! -f "${TRANSCRIPTOME_FASTA}" ]; then
    echo "Generating transcriptome FASTA with gffread ..."
    gffread "${GTF_FILE}" \
        -g "${GENOME_FASTA}" \
        -w "${TRANSCRIPTOME_FASTA}"
    echo "Transcriptome FASTA written to: ${TRANSCRIPTOME_FASTA}"
    echo ""
else
    echo "Transcriptome FASTA already exists at ${TRANSCRIPTOME_FASTA}, skipping"
    echo ""
fi


# --- collect sample prefixes from STAR BAM outputs ---
prefixes=( $( ls ${indir}/bams/*/*.bam | cut -d'A' -f1 ) )
prefixes_uniq=($(echo "${prefixes[@]}" | tr ' ' '\n' | sort -u))

echo "Found ${#prefixes_uniq[@]} samples"
echo ""

# --- run Salmon for each sample ---
for prefix in "${prefixes_uniq[@]}"; do
    sample_name=$( basename "${prefix}" )
    transcriptome_bam="${prefix}Aligned.toTranscriptome.out.bam"

    echo "* Processing sample: ${sample_name}"

    # Check that the transcriptome BAM exists for this sample
    if [ ! -f "${transcriptome_bam}" ]; then
        echo "WARNING: missing transcriptome BAM for ${sample_name}, skipping: ${transcriptome_bam}"
        echo ""
        continue
    fi

    # Skip if Salmon output already exists for this sample
    if [ -f "${SALMON_OUTDIR}/${sample_name}/quant.sf" ]; then
        echo "  quant.sf already exists for ${sample_name}, skipping"
        echo ""
        continue
    fi

    echo "  Input BAM:  ${transcriptome_bam}"
    echo "  Output dir: ${SALMON_OUTDIR}/${sample_name}"
    echo "  RUN"

    # --libType A        : auto-detect library strandedness
    # --alignments       : alignment-based mode using STAR transcriptome BAM
    # --gcBias           : correct for GC bias
    # --seqBias          : correct for sequence-specific bias
    salmon quant \
        -t "${TRANSCRIPTOME_FASTA}" \
        --libType A \
        --alignments "${transcriptome_bam}" \
        --output "${SALMON_OUTDIR}/${sample_name}" \
        --threads "${THREADS}" \
        --gcBias \
        --seqBias    
    echo ""
done

echo "Salmon quantification complete."
echo "Results written to: ${SALMON_OUTDIR}"
echo ""
echo "To import into R with tximport:"
echo "  files <- file.path('${SALMON_OUTDIR}', sample_names, 'quant.sf')"
echo "  txi <- tximport(files, type='salmon', tx2gene=tx2gene_df)"