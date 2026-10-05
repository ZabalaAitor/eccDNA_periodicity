#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-dyad}"
mkdir -p "$OUT_DIR"

GTF="${OUT_DIR}/gencode.v50.annotation.gtf.gz"

wget -q -O "$GTF" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz

gunzip -c "$GTF" |
awk 'BEGIN{OFS="\t"} $3=="gene" {print $1, $4-1, $5}' |
sort -k1,1V -k2,2n > "${OUT_DIR}/hg38_genic_regions.bed"