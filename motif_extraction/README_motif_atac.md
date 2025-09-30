# README: motif_atac_intersect.sh

## Description
This script converts motif bigBed files (`.bb`) into BED format, filters motifs by score, 
intersects them with ATAC peak regions, and generates a new bigBed file for visualization 
(e.g. in GUANACO or IGV).

## Requirements
The following tools must be installed and available in your PATH:

- **UCSC Kent Utilities**  
  * `bigBedToBed`  
  * `bedToBigBed`  
  Download: [UCSC Utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)

- **bedtools** (for intersection)  
  Install with conda:  
  ```bash
  conda install -c bioconda bedtools
  ```

- **Standard Unix utilities** (usually pre-installed on Linux/macOS):  
  `awk`, `sort`, `mkdir`, `mv`

## Inputs
- `motif.bb` : Motif bigBed file (e.g. `JASPAR2024_hg38.bb`)  
- `atac_peaks.bed` : ATAC peak regions in BED format  
- `chrom.sizes` : UCSC chromosome sizes file for the reference genome  
- `fields.as` : Annotation schema file (if required by your bigBed format)  
- `score_threshold` : Minimum motif score for filtering  

## Usage
Example command:

```bash
./motif_atac_intersect.sh JASPAR2024_hg38.bb atac_peaks.bed hg38.chrom.sizes fields.as 400
```

This will output a new bigBed file (`motifs_in_peaks.bb`) containing motifs that intersect 
with ATAC peaks and pass the score threshold.

## Output
- `output/motifs_in_peaks.bb` : Final bigBed file ready for visualization
