# 🧬 Extrachromosomal circular DNA formation follows a 10.5-bp periodicity reflecting the minor groove orientation around nucleosomes

Code and analysis notebook for the study **"Extrachromosomal circular DNA
formation follows a 10.5-bp periodicity reflecting the minor groove orientation
around nucleosomes"**.

The analysis measures the length distribution of eccDNA across 25 published
datasets from four species, and maps eccDNA breakpoints relative to nucleosome
dyads in the four human datasets.

## 📁 Repository structure

```
eccDNA_periodicity.ipynb     Main notebook: runs the pipeline and makes the figures
src/                         Analysis modules and pipeline scripts
environment.yml              Conda environment (Python packages + command-line tools)
requirements.txt             Python dependencies only
data/                        Inputs (downloaded from Zenodo, see below)
results/                     Everything the pipeline produces
```

`src/` contains the Python modules imported by the notebook:

| File | Purpose |
|---|---|
| `preprocessing.py` | Load and merge the BED datasets, dataset-size summaries |
| `length_distribution.py` | Length distributions and KDE peak detection |
| `STFT.py` | Short-time Fourier analysis of the length distributions |
| `dyads.py` | Breakpoint-to-dyad mapping and observed/expected enrichment |
| `call_dyads.py` | Consensus dyad calling from the MNase-seq signal |

and the shell scripts that the notebook calls for the heavy steps:

| Script | Purpose |
|---|---|
| `make_merged_hg38.sh` | Lift the three human Circle-Seq datasets hg19 → GRCh38 and pool them with the microDNA dataset |
| `run_gaffney_consensus_dyads_hg38.sh` | Call consensus nucleosome dyads from the MNase-seq signal and lift them hg18 → GRCh38 |
| `liftover_consensus_dyads_hg18_to_hg38.sh` | LiftOver helper used by the script above |
| `getgenes.sh` | Download GENCODE v50 and derive the genic intervals |
| `non_genic.sh` | Keep only circles that do not overlap a gene |
| `random.sh` | Generate the seeded random permutations |

## ⚙️ Installation

```bash
git clone https://github.com/ZabalaAitor/eccDNA_periodicity
cd eccDNA_periodicity
conda env create -f environment.yml
conda activate eccdna-periodicity
```

The pipeline needs three command-line tools in addition to the Python packages:
`bedtools`, and UCSC `liftOver` and `wigToBigWig`. `environment.yml` installs
all of them; `requirements.txt` covers only the Python side.

> **The `bedtools` version matters.** The pseudo-random stream behind
> `bedtools shuffle -seed` is not stable across versions, so a version other
> than the pinned 2.31.1 produces different permutations and a different
> expected distribution.

## 📥 Input data

Download the input deposit from Zenodo and extract it at the root of the
repository, so that the tree becomes `data/`:

🔗 https://doi.org/10.5281/zenodo.23472112

```
data/
├── arabidopsis/  chicken/  human/  mice/  mice1L3KO/   25 eccDNA datasets (BED3)
└── reference/                                          5 external reference files
```

`data/reference/` holds the MNase-seq midpoint signal (GEO GSE36979), the two
UCSC chain files (hg18→hg38 and hg19→hg38), the GRCh38 chromosome sizes and the
ENCODE blacklist. The GENCODE v50 annotation is **not** part of the deposit:
`src/getgenes.sh` downloads it on demand.

One file in the results deposit is also needed as an input:
`results/dyads/merged_hg38.bed`. The notebook re-derives it and checks the
result, but does not overwrite it, because its row order determines the
`bedtools shuffle -seed` draws and is not reproducible by sorting.

## ▶️ Running the analysis

Open `eccDNA_periodicity.ipynb` and run the cells in order. The pipeline is:

1. pool the human datasets in GRCh38 (`make_merged_hg38.sh`);
2. call consensus nucleosome dyads (`run_gaffney_consensus_dyads_hg38.sh`);
3. derive the genic intervals and keep the non-genic circles
   (`getgenes.sh`, `non_genic.sh`);
4. generate 1000 seeded permutations (`random.sh`);
5. map breakpoints to the closest dyad and compute the observed/expected
   enrichment (`dyads.py`).

Steps 2 and 5 are the slow ones; the dyad calling takes about ten minutes and
the permutation analysis considerably longer.

## 📊 Results

The full set of outputs is deposited at Zenodo:

🔗 https://doi.org/10.5281/zenodo.23474862

## 📄 License

MIT. See [LICENSE](LICENSE).
```
