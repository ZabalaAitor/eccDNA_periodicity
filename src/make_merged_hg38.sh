#!/usr/bin/env bash
set -euo pipefail

# Rebuild the pooled human eccDNA coordinates in hg38 and verify them against the
# deposited copy.
# Usage: bash src/make_merged_hg38.sh [human_dir] [chain_file] [out_dir]
#
# merged_hg38.bed pools the four human datasets in hg38 coordinates:
#   Hs-cancer_cells-MD                      called against GRCh38 by nf-core/circdna, used as is
#   Hs-blood-CS, Hs-muscle-CS, Hs-sperm-CS  published in hg19, lifted over (76 of 22,398 lost)
#
# The regenerated file is written next to the deposited one and only compared with
# it. The row order of the deposited copy fixes the per-interval draws of
# `bedtools shuffle -seed`, so rewriting it in place would change the permutations
# and therefore the enrichment statistics.

HUMAN_DIR="${1:-data/human}"
CHAIN_FILE="${2:-data/reference/hg19ToHg38.over.chain.gz}"
OUT_DIR="${3:-results/dyads}"

LIFTOVER_DIR="${OUT_DIR}/liftover"
REFERENCE_BED="${OUT_DIR}/merged_hg38.bed"
REGENERATED_BED="${OUT_DIR}/merged_hg38.regenerated.bed"
HG19_CHR1_SIZE=249250621
LIFTED_DATASETS=(Hs-blood-CS Hs-muscle-CS Hs-sperm-CS)
NATIVE_DATASET=Hs-cancer_cells-MD

if ! command -v liftOver >/dev/null 2>&1; then
    echo "Error: liftOver not found in PATH." >&2
    exit 1
fi

if [[ ! -f "$CHAIN_FILE" ]]; then
    echo "Error: chain file not found: $CHAIN_FILE" >&2
    exit 1
fi

# Guard: the chain must be hg19 -> hg38. liftOver accepts any chain without
# checking the assembly, so a wrong one fails silently; the chr1 target size is
# unambiguous (hg19 249250621, hg18 247249719, hg38 248956422).
# No early `exit` in awk: it would SIGPIPE zcat and trip `set -o pipefail`.
CHR1_SIZE=$(zcat "$CHAIN_FILE" | awk '$1=="chain" && $3=="chr1" && !seen {print $4; seen=1}')
if [[ "$CHR1_SIZE" != "$HG19_CHR1_SIZE" ]]; then
    echo "Error: $CHAIN_FILE is not an hg19 chain (chr1 target size $CHR1_SIZE)." >&2
    exit 1
fi
echo "Chain OK: hg19 -> hg38"

mkdir -p "$LIFTOVER_DIR"
for NAME in "${LIFTED_DATASETS[@]}"; do
    liftOver \
        "${HUMAN_DIR}/${NAME}.bed" \
        "$CHAIN_FILE" \
        "${LIFTOVER_DIR}/${NAME}_hg38.bed" \
        "${LIFTOVER_DIR}/${NAME}_unmapped.bed"
done
wc -l "${LIFTOVER_DIR}"/*_hg38.bed

cat "${LIFTOVER_DIR}"/Hs-{blood,muscle,sperm}-CS_hg38.bed \
    "${HUMAN_DIR}/${NATIVE_DATASET}.bed" \
    | cut -f1-3 \
    | sort -k1,1 -k2,2n \
    > "$REGENERATED_BED"

if [[ ! -f "$REFERENCE_BED" ]]; then
    echo "Error: $REFERENCE_BED not found." >&2
    echo "It is a provided input of the Zenodo deposit and cannot be substituted by" >&2
    echo "$REGENERATED_BED, whose row order differs and would change the permutations." >&2
    exit 1
fi

EXPECTED_MD5=$(sort -k1,1 -k2,2n "$REFERENCE_BED" | md5sum | cut -d" " -f1)
OBTAINED_MD5=$(md5sum "$REGENERATED_BED" | cut -d" " -f1)
if [[ "$EXPECTED_MD5" != "$OBTAINED_MD5" ]]; then
    echo "Error: the liftOver output does not reproduce $REFERENCE_BED" >&2
    echo "  deposited (sorted) md5: $EXPECTED_MD5" >&2
    echo "  regenerated        md5: $OBTAINED_MD5" >&2
    exit 1
fi

echo "OK: $REFERENCE_BED reproduced from the raw BED files"
echo "    content md5 $OBTAINED_MD5, $(wc -l < "$REFERENCE_BED" | tr -d ' ') intervals"
