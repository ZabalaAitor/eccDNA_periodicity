#!/usr/bin/env bash
set -euo pipefail

# Call consensus nucleosome dyads from MNase-seq midpoint signal and lift them to hg38.
# Usage: bash src/run_gaffney_consensus_dyads_hg38.sh <mids.wig.gz> <hg18ToHg38.chain.gz> <out.tsv>
# Converts the WIG track to BigWig, scores positions with call_dyads.py, then lifts hg18 -> hg38.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

INPUT_WIG_GZ="${1:-${PROJECT_DIR}/dyad/GSE36979_mnase_mids_combined_147.wig.gz}"
CHAIN_FILE="${2:-${PROJECT_DIR}/dyad/hg18ToHg38.over.chain.gz}"
OUTPUT_TSV="${3:-${PROJECT_DIR}/dyad/gaffney_consensus_dyads_hg38.tsv}"

# Convert input paths to absolute paths before changing directory.
INPUT_WIG_GZ="$(realpath "$INPUT_WIG_GZ")"
CHAIN_FILE="$(realpath "$CHAIN_FILE")"
OUTPUT_TSV="$(realpath "$OUTPUT_TSV")"

if [[ ! -f "$INPUT_WIG_GZ" ]]; then
    echo "Error: input WIG.GZ not found: $INPUT_WIG_GZ" >&2
    exit 1
fi

if [[ ! -f "$CHAIN_FILE" ]]; then
    echo "Error: chain file not found: $CHAIN_FILE" >&2
    exit 1
fi

if ! command -v wigToBigWig >/dev/null 2>&1; then
    echo "Error: wigToBigWig not found in PATH." >&2
    exit 1
fi

if ! command -v liftOver >/dev/null 2>&1; then
    echo "Error: liftOver not found in PATH." >&2
    exit 1
fi

TMP_DIR="$(mktemp -d)"

cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# The hg18 chromosome sizes are read from the chain file itself: every `chain`
# header line carries the target (hg18) chromosome name and its full length.
# This keeps the step self-contained -- no download is needed -- and guarantees
# that the sizes always match the chain actually in use.
CHROM_SIZES_TMP="${TMP_DIR}/hg18.chrom.sizes"

gunzip -c "$CHAIN_FILE" \
    | awk '$1 == "chain" { print $3 "\t" $4 }' \
    | sort -u \
    > "$CHROM_SIZES_TMP"

if [[ ! -s "$CHROM_SIZES_TMP" ]]; then
    echo "Error: no hg18 chromosome sizes could be read from $CHAIN_FILE." >&2
    exit 1
fi

cp "$INPUT_WIG_GZ" \
    "$TMP_DIR/GSE36979_mnase_mids_combined_147.wig.gz"

gunzip -c "$TMP_DIR/GSE36979_mnase_mids_combined_147.wig.gz" \
    > "$TMP_DIR/GSE36979_mnase_mids_combined_147.wig"

wigToBigWig \
    "$TMP_DIR/GSE36979_mnase_mids_combined_147.wig" \
    "$CHROM_SIZES_TMP" \
    "$TMP_DIR/GSE36979_mnase_mids_combined_147.bw"

(
    cd "$TMP_DIR"
    python3 "$SCRIPT_DIR/call_dyads.py"
)

(
    cd "$TMP_DIR"
    bash "$SCRIPT_DIR/liftover_consensus_dyads_hg18_to_hg38.sh" \
        "gaffney_consensus_dyads_hg18.tsv" \
        "$CHAIN_FILE" \
        "gaffney_consensus_dyads_hg38"
)

mv \
    "$TMP_DIR/gaffney_consensus_dyads_hg38.tsv" \
    "$OUTPUT_TSV"

echo "Done. Final output: $OUTPUT_TSV"