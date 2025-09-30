#!/bin/bash

# ================================
# Script: Motif–ATAC Peak Intersection
# Purpose: Convert motif bigBed (.bb) into bed, filter by score,
#          intersect with ATAC peaks, and generate a new bigBed for visualization.
# Requirements: UCSC tools (bigBedToBed, bedToBigBed), bedtools, genome chrom.sizes file
# ================================

# ---- Usage ----
# ./motif_atac_intersect.sh <motif.bb> <atac_peaks.bed> <chrom.sizes> <fields.as> <score_threshold>
#
# Example:
# ./motif_atac_intersect.sh JASPAR2024_hg38.bb atac_peaks.bed hg38.chrom.sizes fields.as 400
#
# Output: output/motifs_in_peaks.bb

# ---- Input arguments ----
MOTIF_BB=$1          # motif bigBed file
ATAC_PEAK_BED=$2     # ATAC peak regions in BED format
CHROM_SIZES=$3       # UCSC chrom.sizes file for your genome
FIELDS_AS=$4         # annotation schema (if required)
THRESHOLD=$5         # motif score threshold

# ---- Check inputs ----
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <motif.bb> <atac_peaks.bed> <chrom.sizes> <fields.as> <score_threshold>"
    exit 1
fi

# ---- Intermediate and output names ----
RAW_BED="motifs_raw.bed"
FILTERED_BED="motifs_filtered.bed"
INTERSECT_BED="motifs_atac_intersect.bed"
INTERSECT_SORTED="motifs_atac_intersect_sorted.bed"
OUTPUT_BB="motifs_in_peaks.bb"
OUTPUT_DIR="output"

# ================================
# Step 1: Convert bigBed to bed
# ================================
echo "[1/4] Converting motif .bb to .bed ..."
bigBedToBed "$MOTIF_BB" "$RAW_BED"

# ================================
# Step 2: Filter motifs by score
# ================================
echo "[2/4] Filtering motifs (score > $THRESHOLD) ..."
awk -v thr="$THRESHOLD" '$5 > thr' "$RAW_BED" > "$FILTERED_BED"

# ================================
# Step 3: Intersect filtered motifs with ATAC peaks
# ================================
echo "[3/4] Intersecting motifs with ATAC peaks ..."
bedtools intersect -a "$FILTERED_BED" -b "$ATAC_PEAK_BED" > "$INTERSECT_BED"

# Sort the intersections
echo "Sorting intersections ..."
sort -k1,1 -k2,2n "$INTERSECT_BED" > "$INTERSECT_SORTED"

# ================================
# Step 4: Convert intersected bed to bigBed
# ================================
echo "[4/4] Converting intersected motifs to .bb ..."
mkdir -p "$OUTPUT_DIR"
bedToBigBed -as="$FIELDS_AS" -type=bed6+1 "$INTERSECT_SORTED" "$CHROM_SIZES" "$OUTPUT_BB"
mv "$OUTPUT_BB" "$OUTPUT_DIR/"

echo "✅ Done. Output written to: $OUTPUT_DIR/$OUTPUT_BB"

