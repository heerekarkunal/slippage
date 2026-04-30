# PCR Slippage Artefact Detection Pipeline

Detects, locates, and classifies insertion/deletion mutations in Illumina MiSeq
NGS data as either **PCR-induced slippage artefacts** or **genuine biological
mutations**, using reference-guided homopolymer/STR mapping and Phred quality
scoring.

---

## What this pipeline does

```
Raw MiSeq data                    Reference FASTA
(BCL or FASTQ)                         │
       │                               │
       ▼                               ▼
  [BCL→FASTQ]             [Build slippage hotspot map]
       │                    (homopolymers + STRs)
       │                               │
       ▼                               │
  [Align reads]                        │
  bwa-mem2 → BAM                       │
       │                               │
       ▼                               │
  [Call indels]                        │
  bcftools mpileup/call                │
       │                               │
       └──────────┬────────────────────┘
                  ▼
         [Classify each indel]
                  │
     ┌────────────┼──────────────┐
     ▼            ▼              ▼
POTENTIAL   HIGHLY POTENTIAL  LOW POTENTIAL
PCR         NATURAL          NATURAL
ARTEFACT    MUTATION         MUTATION
     │
     └── overlaps homopolymer or STR in reference
```

### Classification logic

| Label | Condition |
|---|---|
| **POTENTIAL PCR ARTEFACT** | Indel falls within a homopolymer or STR region of the reference |
| **HIGHLY POTENTIAL NATURAL MUTATION** | Indel outside any hotspot AND mean supporting-read Phred ≥ threshold |
| **LOW POTENTIAL NATURAL MUTATION** | Indel outside any hotspot AND mean supporting-read Phred < threshold |

---

## Installation

### 1. Clone/download the project

```bash
git clone <your-repo-url>
cd slippage_pipeline
```

### 2. Set up a conda environment (recommended)

```bash
conda create -n slippage python=3.11
conda activate slippage
```

### 3. Install system bioinformatics tools

```bash
# All required tools via conda (easiest, one command):
conda install -c bioconda -c conda-forge \
    bwa-mem2 samtools bcftools tabix htslib

# Or via apt (Ubuntu/Debian):
sudo apt update
sudo apt install bwa samtools bcftools tabix
# Note: apt usually has older versions. conda is preferred.
```

For BCL input mode only — install bcl-convert:
```
Download from: https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html
```

### 4. Install Python packages

```bash
pip install -r requirements.txt
```

### 5. Verify everything is ready

```bash
# Quick dependency check (the pipeline runs this automatically,
# but you can check manually):
samtools --version
bcftools --version
bwa-mem2 version
python -c "import pysam, biopython, pandas, numpy, matplotlib, regex; print('OK')"
```

---

## Usage

### Basic usage — from FASTQ

```bash
python run_pipeline.py \
    --fastq sample.fastq.gz \
    --reference genome.fasta \
    --output ./results
```

### Paired-end FASTQ

```bash
python run_pipeline.py \
    --fastq1 sample_R1.fastq.gz \
    --fastq2 sample_R2.fastq.gz \
    --reference genome.fasta \
    --output ./results
```

### From raw BCL run folder

```bash
python run_pipeline.py \
    --bcl /path/to/MiSeqRunFolder \
    --reference genome.fasta \
    --output ./results
```

### Full example with custom thresholds

```bash
python run_pipeline.py \
    --fastq1 reads_R1.fastq.gz \
    --fastq2 reads_R2.fastq.gz \
    --reference my_organism.fasta \
    --output ./results \
    --min-hp-len 6 \
    --min-str-copies 3 \
    --min-quality 20 \
    --high-quality-threshold 30 \
    --min-depth 10 \
    --min-vaf 0.02 \
    --threads 8 \
    --keep-intermediates \
    --verbose
```

---

## All options

```
Input (choose one mode):
  --fastq FILE          Single-end FASTQ (.fastq or .fastq.gz)
  --fastq1 FILE         Paired-end Read 1
  --fastq2 FILE         Paired-end Read 2
  --bcl DIR             Raw MiSeq BCL run folder

Required:
  --reference FASTA     Reference genome FASTA
  --output DIR          Output directory

Slippage hotspot thresholds:
  --min-hp-len N        Minimum homopolymer length to flag (default: 6)
  --min-str-unit N      Minimum STR unit length in bp (default: 2)
  --max-str-unit N      Maximum STR unit length in bp (default: 6)
  --min-str-copies N    Minimum STR copy number to flag (default: 3)

Mutation classification:
  --min-quality Q           Minimum mean read Phred score (default: 20)
  --high-quality-threshold  Phred score for HIGHLY POTENTIAL NATURAL cutoff (default: 30)
  --min-depth N             Minimum read depth at variant site (default: 5)
  --min-vaf F               Minimum variant allele frequency 0–1 (default: 0.01)

Alignment:
  --aligner             bwa-mem2 / bwa / minimap2 (default: bwa-mem2)
  --threads N           CPU threads (default: 4)

BCL:
  --sample-sheet FILE   Override SampleSheet.csv path

Misc:
  --keep-intermediates  Retain BAM and intermediate VCF files
  --verbose             Verbose console output
```

---

## Output files

After a successful run, `results/reports/` contains:

```
reports/
├── variants_classified.tsv      ← Full table of every indel, annotated
├── per_site_counts.tsv          ← Count of insertions + deletions at each position
├── artefacts_flagged.tsv        ← Only PCR artefact calls
├── natural_mutations.tsv        ← Only natural mutation calls
├── hotspot_map_summary.tsv      ← All slippage-prone regions in the reference
├── summary_report.txt           ← Human-readable text summary
└── figures/
    ├── classification_overview.png   ← Pie + bar chart of classifications
    ├── quality_distribution.png      ← Violin plot of quality per class
    ├── genome_map.png                ← Chromosome-scale indel position map
    └── hotspot_positions.png         ← Hotspot density across the genome
```

If `--keep-intermediates` is used, `results/intermediates/` also contains:
```
intermediates/
├── *.sorted.bam         ← Aligned, sorted reads
├── *.sorted.bam.bai     ← BAM index
├── *.indels.vcf.gz      ← Called, normalised indels
├── *.indels.vcf.gz.tbi  ← VCF index
├── hotspots.bed         ← Hotspot regions (load in IGV)
├── hotspots.tsv         ← Hotspot regions (human-readable)
└── *.flagstat.txt       ← Alignment statistics
```

### Key output columns in `variants_classified.tsv`

| Column | Description |
|---|---|
| `chrom` | Contig/chromosome name |
| `pos_1` | 1-based genomic position |
| `ref` / `alt` | Reference and alternate alleles |
| `mut_type` | `insertion` or `deletion` |
| `mut_size_bp` | Size of the indel in bp |
| `total_depth` | Total read depth at this site |
| `vaf_percent` | Variant allele frequency (%) |
| `mean_alt_supporting_phred` | Mean Phred quality of ALT-supporting reads |
| `in_hotspot` | True/False |
| `hotspot_type` | `homopolymer` or `STR` |
| `hotspot_unit` | Repeat unit (e.g., `A` or `AT`) |
| `hotspot_risk_level` | `HIGH` or `MODERATE` |
| `classification` | The final label |
| `classification_reason` | Detailed explanation of why this label was assigned |

---

## Visualising results in IGV

The `hotspots.bed` file in `intermediates/` can be loaded directly into
[IGV (Integrative Genomics Viewer)](https://igv.org/):

1. Open IGV and load your reference FASTA
2. Load the sorted BAM file
3. Load `hotspots.bed` — hotspot regions will appear as coloured tracks
4. Load the VCF file to see variant calls overlaid on coverage

---

## Adjusting thresholds for your organism

For a ~7300 bp organism (viral/bacterial/plasmid), these defaults are appropriate:

- `--min-hp-len 6` : catches all biologically relevant homopolymers. Increase to 8
  if you want only the highest-risk runs (reduces false positive artefact calls).
- `--min-str-copies 3` : 3 copies of any 2–6 bp unit. Sufficient for a small genome.
- `--high-quality-threshold 30` : Q30 = 99.9% base call accuracy. The standard
  cutoff in clinical NGS. Use Q25 if your run has generally lower quality.
- `--min-vaf 0.01` : 1% VAF. Lower this (e.g. 0.005) if looking for rare variants;
  raise it (e.g. 0.05) to focus on high-frequency variants only.

---

## Interpreting results

**If most of your indels are `POTENTIAL PCR ARTEFACT`:**
Your amplicon likely passes through repetitive sequence. This is expected for
many organisms. Consider: (1) validating key sites by Sanger sequencing,
(2) using a higher-fidelity polymerase (Q5, Phusion) to reduce stutter,
or (3) increasing PCR cycle numbers cautiously.

**If you have `HIGHLY POTENTIAL NATURAL MUTATION` calls:**
These are your best candidates for genuine biological mutations. Check the
`vaf_percent` column — high-VAF calls (>30%) are very likely real; low-VAF
calls (1–5%) should be confirmed by an independent method.

**If most calls are `LOW POTENTIAL NATURAL MUTATION`:**
The run had low overall quality. Reprocess with stricter quality filters
or re-sequence with higher coverage.
