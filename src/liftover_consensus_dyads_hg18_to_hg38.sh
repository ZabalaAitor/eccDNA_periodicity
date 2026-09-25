#!/usr/bin/env bash
set -euo pipefail

# Lift over consensus dyads from hg18 to hg38.
# Input TSV columns expected:
# chrom start end <...remaining columns...>
# Only columns 1-3 (chrom/start/end) are lifted; all other columns are appended.

INPUT_TSV="${1:-gaffney_consensus_dyads_hg18.tsv}"
CHAIN_FILE="${2:-hg18ToHg38.over.chain.gz}"
OUTPUT_PREFIX="${3:-gaffney_consensus_dyads_hg38}"

TMP_BED="${OUTPUT_PREFIX}.tmp.hg18.bed"
TMP_ROWS="${OUTPUT_PREFIX}.tmp.rows.tsv"
MAPPED_BED="${OUTPUT_PREFIX}.bed"
UNMAPPED_BED="${OUTPUT_PREFIX}.unmapped.bed"
MAPPED_TSV="${OUTPUT_PREFIX}.tsv"
UNMAPPED_TSV="${OUTPUT_PREFIX}.unmapped.tsv"

if ! command -v liftOver >/dev/null 2>&1; then
    echo "Error: liftOver not found in PATH." >&2
    exit 1
fi

if [[ ! -f "$INPUT_TSV" ]]; then
    echo "Error: input file not found: $INPUT_TSV" >&2
    exit 1
fi

if [[ ! -f "$CHAIN_FILE" ]]; then
    echo "Error: chain file not found: $CHAIN_FILE" >&2
    exit 1
fi

# Convert TSV to strict BED6 for liftOver compatibility.
# Only columns 1-3 (chrom/start/end) are used for liftover.
# All other columns are preserved by row ID and restored afterward.
# BED6: chrom, start, end, name, score, strand
# Note: liftOver can misinterpret files with >=12 columns as BED12.
awk 'BEGIN { OFS="\t" }
NR == 1 { next }
{
    if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) {
        printf("Error: non-integer start/end at input line %d: %s\n", NR, $0) > "/dev/stderr"
        exit 2
    }
    if ($2 >= $3) {
        printf("Error: start >= end at input line %d: %s\n", NR, $0) > "/dev/stderr"
        exit 2
    }
    id = "row_" NR
    print $1, $2, $3, id, 0, "."
}' "$INPUT_TSV" > "$TMP_BED"

# Write the ID->row mapping used to restore full TSV columns.
awk 'BEGIN { OFS="\t" }
NR == 1 { next }
{
    id = "row_" NR
    print id, $0
}' "$INPUT_TSV" > "$TMP_ROWS"

# Run liftOver.
liftOver "$TMP_BED" "$CHAIN_FILE" "$MAPPED_BED" "$UNMAPPED_BED"

HEADER="$(head -n 1 "$INPUT_TSV")"

# Rebuild mapped TSV with lifted coordinates and all original columns.
{
    printf '%s\n' "$HEADER"
    awk 'BEGIN { OFS="\t" }
    FNR == NR {
        id = $1
        $1 = ""
        sub(/^\t/, "")
        row[id] = $0
        next
    }
    {
        id = $4
        if (!(id in row)) {
            next
        }
        n = split(row[id], fields, "\t")
        fields[1] = $1
        fields[2] = $2
        fields[3] = $3
        out = fields[1]
        for (i = 2; i <= n; i++) {
            out = out OFS fields[i]
        }
        print out
    }' "$TMP_ROWS" "$MAPPED_BED"
} > "$MAPPED_TSV"

# Rebuild unmapped TSV with original coordinates and all original columns.
{
    printf '%s\n' "$HEADER"
    awk 'BEGIN { OFS="\t" }
    FNR == NR {
        id = $1
        $1 = ""
        sub(/^\t/, "")
        row[id] = $0
        next
    }
    /^#/ { next }
    {
        id = $4
        if (id in row) {
            print row[id]
        }
    }' "$TMP_ROWS" "$UNMAPPED_BED"
} > "$UNMAPPED_TSV"

# Clean temp files.
rm -f "$TMP_BED" "$TMP_ROWS"

echo "LiftOver complete."
echo "Mapped BED:                 $MAPPED_BED"
echo "Unmapped BED:               $UNMAPPED_BED"
echo "Mapped TSV (all columns):   $MAPPED_TSV"
echo "Unmapped TSV (all columns): $UNMAPPED_TSV"
