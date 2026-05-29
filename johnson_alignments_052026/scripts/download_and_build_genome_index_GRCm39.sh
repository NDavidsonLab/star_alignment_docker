#!/bin/bash

script_dir=$(dirname "$(realpath "$0")")
cd "${script_dir}/.."

GENOMEDIR="./genome/M37"

mkdir -p ${GENOMEDIR}

(
    cd ${GENOMEDIR}
    # download the "Genome sequence, primary assembly (GRCm39)" fasta file
    [ -f GRCm39.primary_assembly.genome.fa.gz ] || \
        wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/GRCm39.primary_assembly.genome.fa.gz
    # filter it as described in the note below
    # download the annotations that correspond to it 
    [ -f gencode.vM37.primary_assembly.annotation.gtf.gz ] || \
        wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.primary_assembly.annotation.gtf.gz

    gunzip -c GRCm39.primary_assembly.genome.fa.gz > GRCm39.primary_assembly.genome.fa
    gunzip -c gencode.vM37.primary_assembly.annotation.gtf.gz > gencode.vM37.primary_assembly.annotation.gtf

)

mkdir -p ${GENOMEDIR}/STAR
STAR --runThreadN 4 --runMode genomeGenerate \
    --genomeDir ${GENOMEDIR}/STAR \
    --genomeFastaFiles ${GENOMEDIR}/GRCm39.primary_assembly.genome.fa \
    --sjdbGTFfile ${GENOMEDIR}/gencode.vM37.primary_assembly.annotation.gtf.gz \
