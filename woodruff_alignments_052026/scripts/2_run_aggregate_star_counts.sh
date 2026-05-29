#!/bin/bash

# --- configuration ---

script_dir=$(dirname "$(realpath "$0")")
cd "${script_dir}/../"


# inputs
genome_dir="./genome/M37/STAR"
data_dir="./data"
indir="./alignments"

outdir="./abundances"
mkdir -p ${outdir}
mkdir -p ${outdir}/STAR

COUNT_COL=2       # column to extract: 2=unstranded, 3=stranded forward, 4=stranded reverse
OUTFILE="${outdir}/STAR/counts_matrix.tsv"

# Aggregates STAR ReadsPerGene.out.tab count files across all samples into a
# single tab-separated count matrix. STAR outputs 4 columns per file:
#   col 1: gene ID
#   col 2: unstranded counts
#   col 3: stranded counts (forward)
#   col 4: stranded counts (reverse)
# This script extracts the unstranded counts (col 2) by default -- change the
# column index below if your library is stranded.
#
# Output: ${outdir}/counts_matrix.tsv
#   - rows are genes
#   - columns are samples (named by prefix)
#   - first 4 rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
#     are STAR summary stats and can be removed if needed


# --- collect sample prefixes ---
prefixes=( $( ls ${indir}/bams/*/*.bam | cut -d'A' -f1 ) )
prefixes_uniq=($(echo "${prefixes[@]}" | tr ' ' '\n' | sort -u))

echo "Found ${#prefixes_uniq[@]} samples"

# --- build the matrix ---
# Start with the gene IDs from the first sample as the first column
first_prefix="${prefixes_uniq[0]}"
first_file="${first_prefix}ReadsPerGene.out.tab"

if [ ! -f "${first_file}" ]; then
    echo "ERROR: could not find ${first_file}"
    exit 1
fi

# Extract gene IDs to a tmp file (skip first 4 summary rows)
tmp_matrix=$(mktemp)
awk 'NR>4 {print $1}' "${first_file}" > "${tmp_matrix}"

# Build a header line starting with "ENSG_ID"
header="ENSG_ID"

# For each sample, paste its count column alongside the gene IDs
for prefix in "${prefixes_uniq[@]}"; do
    count_file="${prefix}ReadsPerGene.out.tab"
    sample_name=$( basename "${prefix}" )

    if [ ! -f "${count_file}" ]; then
        echo "WARNING: missing count file for ${sample_name}, skipping: ${count_file}"
        continue
    fi

    echo "Adding sample: ${sample_name}"

    # Extract the chosen count column and paste it to the growing matrix
    tmp_counts=$(mktemp)
    awk -v col="${COUNT_COL}" 'NR>4 {print $col}' "${count_file}" > "${tmp_counts}"
    paste "${tmp_matrix}" "${tmp_counts}" > "${tmp_matrix}.new"
    mv "${tmp_matrix}.new" "${tmp_matrix}"
    rm "${tmp_counts}"

    header="${header}\t${sample_name}"
done

# Write the header and matrix to the output file
echo -e "${header}" | cat - "${tmp_matrix}" > "${OUTFILE}"

rm "${tmp_matrix}"

echo ""
echo "Count matrix written to: ${OUTFILE}"
echo "Dimensions: $(wc -l < ${OUTFILE}) rows x $(head -1 ${OUTFILE} | tr '\t' '\n' | wc -l) columns"