#!/usr/bin/env bash
set -euo pipefail

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

CHROM_SIZES_LOCAL="${SCRIPT_DIR}/hg18.chrom.sizes"
CHROM_SIZES_TMP="${TMP_DIR}/hg18.chrom.sizes"

if [[ -f "$CHROM_SIZES_LOCAL" ]]; then
    cp "$CHROM_SIZES_LOCAL" "$CHROM_SIZES_TMP"
else
    if command -v wget >/dev/null 2>&1; then
        wget -q -O "$CHROM_SIZES_TMP" \
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/bigZips/hg18.chrom.sizes"
    elif command -v curl >/dev/null 2>&1; then
        curl -fsSL \
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/bigZips/hg18.chrom.sizes" \
            -o "$CHROM_SIZES_TMP"
    else
        echo "Error: need wget or curl to download hg18.chrom.sizes." >&2
        exit 1
    fi
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