# 🧬 Extrachromosomal circular DNA formation follows a 10.5-bp periodicity imposed by the minor groove orientation around nucleosomes

This repository contains the code and notebooks used in the study **"Extrachromosomal circular DNA formation follows a 10.5 base-pair periodicity imposed by the minor groove orientation around nucleosomes"**.

## 📁 Repository Structure

- `eccDNA_periodicity_paper.ipynb`: Main notebook performing the analysis and generating plots.

- `src`: Script for running the notebook.

## 📥 How to Use This Repository

### 1. Clone the Repository

```bash
git clone https://github.com/ZabalaAitor/eccDNA_periodicity
cd eccDNA_periodicity
```

### 2. Download input data from Zenodo:

🔗 [https://zenodo.org/records/15783792](https://zenodo.org/record/15783792)

After downloading, extract the contents and organize them as follows:

- `data/` → should contain all BED files with the eccDNA detection.

- `dyad/` → contains nucleosome dyad positions and related reference files used for nucleosome positioning/periodicity analysis:
  - `gaffney_consensus_dyads_hg38.tsv` → consensus nucleosome dyad positions in the human hg38 genome assembly.
  - `GSE36979_mnase_mids_combined_147.wig.gz` → MNase-seq nucleosome dyad/midpoint signal, corresponding to nucleosome positioning data.
  - `hg18ToHg38.over.chain.gz` → UCSC chain file used to convert genomic coordinates from hg18 to hg38.

These two folders (`data/` and `dyad/`) are essential for running the pipeline.

## 📊 Outputs

The ouputs of the benchmarking analysis are available at:  
🔗 [https://zenodo.org/records/15783794](https://zenodo.org/records/15783794)

## Cite
