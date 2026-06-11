# Phylogeographic analysis of Staphylococcus nepalensis reveals global occurrence of antimicrobial-resistant lineages carrying the intrinsic sal(E) resistance gene

## Overview

This repository contains the complete bioinformatics pipeline used for the comparative genomic analysis of *Staphylococcus nepalensis*, including genome quality assessment, annotation, pangenome analysis, phylogenetic inference, and antimicrobial resistance (AMR) profiling. The pipeline processes 36 high-quality genomes to identify the global phylogeographic structure of this coagulase-negative staphylococcus species and characterizes the genetic basis of the intrinsic sal(E) lincosamide resistance gene.

**Key findings:**
- All 36 *S. nepalensis* genomes carry the intrinsic sal(E) gene (100% core resistome)
- Single ancestral acquisition with minimal subsequent modification (97.8% sequence identity across 543 amino acids)
- Accessory AMR genes show variable distribution: fosB (58.3%), tet(K) (25.0%), mph(C) (16.7%), mecA (5.6%)
- Plasmid-associated replicons (rep7a, rep7b, rep19c, rep21) link genetic determinants to mobile elements
- No virulence genes detected across all genomes, consistent with non-pathogenic CoNS status

---

## Pipeline Overview

```
1. NCBI Genome Download
   ├─ Tool: NCBI Datasets CLI v18.29.1
   ├─ Query: Taxon "Staphylococcus nepalensis"
   └─> 88 genomes downloaded

2. Deduplication by BioSample
   ├─ Groups genomes by BioSample accession
   ├─ Keeps highest assembly version per BioSample
   ├─ Prefers GCF over GCA accessions
   ├─ 42 duplicates removed
   └─> 46 unique genomes

3. Quality Filtering (CheckM2 v1.1.0 + QUAST v5.2.0)
   ├─ Completeness >= 95%         (CheckM2)
   ├─ Contamination <= 5%         (CheckM2)
   ├─ Total contigs <= 100        (QUAST)
   ├─ Total N's <= 500            (QUAST)
   ├─ Exclude suppressed genomes  (NCBI status)
   ├─ 10 genomes excluded (low quality or suppressed)
   └─> 36 high-quality genomes retained

4. Parallel Processing
   ├─ FastANI v1.34
   │  ├─ All-vs-all ANI matrix (36 × 36)
   │  ├─ All-vs-reference (GCF_002442935.1 = J11 reference)
   │  └─> ANI range: 98.99% – 99.xx% (all > 95% species threshold)
   │
   ├─ Bakta v1.12.0 (full BaktaDB v6.0)
   │  ├─ Annotation of all 36 high-quality genomes
   │  ├─ Output: .gff3, .gbff, .faa, .ffn, .tsv per genome
   │  └─> Used as input for Panaroo and AMRFinderPlus
   │
   ├─ Panaroo v1.3.4
   │  ├─ Input: 36 annotated .gff3 files from Bakta
   │  ├─ Mode: strict clean-mode
   │  ├─ Core threshold: 98% (≥98% of isolates)
   │  ├─ Core genes:       2,348
   │  ├─ Soft core genes:     40
   │  ├─ Shell genes:         485
   │  ├─ Cloud genes:       1,880
   │  ├─ Total pangenome:   4,753
   │  ├─ Pangenome type:    Open
   │  └─> core_gene_alignment_filtered.aln → input for Gubbins
   │
   └─ VirulenceFinder + Abricate v(VFDB)
      ├─ Databases: s.aureus_toxin, s.aureus_exoenzyme,
      │             s.aureus_hostimm, VFDB
      ├─ Identity threshold: ≥90% (VirulenceFinder), ≥80% (Abricate)
      ├─ Coverage threshold: ≥60%
      └─> No virulence genes detected in any of 36 genomes

5. AMR Analysis
   ├─ AMRFinderPlus v4.2.7 (DB: 2026-05-15.1)
   │  ├─ Organism: Staphylococcus aureus
   │  ├─ Flag: --plus (enables stress and heavy metal gene detection)
   │  ├─ Input: .faa protein files from Bakta
   │  │
   │  ├─ Core resistome (100% genomes):
   │  │  ├─ sal(E)  -- Lincosamide/Pleuromutilin/Streptogramin A
   │  │
   │  ├─ Accessory AMR genes:
   │  │  ├─ fosB/fosB4          21/36 (58.3%) -- Fosfomycin
   │  │  ├─ tet(K)               9/36 (25.0%) -- Tetracycline
   │  │  ├─ mph(C)               6/36 (16.7%) -- Macrolide
   │  │  ├─ qacG                 5/36 (13.9%) -- Quaternary ammonium (biocide)
   │  │  ├─ str                  4/36 (11.1%) -- Aminoglycoside
   │  │  ├─ lnu(A)/lnu(A)'       3/36  (8.3%) -- Lincosamide
   │  │  ├─ catA                 3/36  (8.3%) -- Phenicol
   │  │  ├─ aac(6')-Ie/aph(2'') 2/36  (5.6%) -- Aminoglycoside
   │  │  ├─ mecA/mecI/mecR1      2/36  (5.6%) -- Beta-lactam/MRSA
   │  │  ├─ dfrE                 2/36  (5.6%) -- Trimethoprim
   │  │  ├─ erm(B)               2/36  (5.6%) -- Macrolide
   │  │  ├─ dfrG                 1/36  (2.8%) -- Trimethoprim
   │  │  └─ blaI/blaPC1/blaR1    1/36  (2.8%) -- Beta-lactam regulators
   │  │
   │  └─ Heavy metal & stress resistance genes (via --plus flag):
   │     ├─ arsR/arsB/arsC  -- Arsenic resistance operon (widely distributed)
   │     ├─ cadD             -- Cadmium tolerance
   │     ├─ mco              -- Copper tolerance / oxidative stress
   │     └─ merA/merB/merT   -- Mercury resistance operon
   │
   ├─ RGI/CARD (DB v3.2.7)
   │  ├─ Input: .faa protein files from Bakta
   │  ├─ Reporting: Perfect and Strict hits only
   │  └─> Confirmed: sal(E), tet(K), mecA, lnu(A), mph(C), qacG
   │
   ├─ PlasmidFinder v2.2.0 (DB: 2024-01-24)
   │  ├─ Input: .fna genome assemblies
   │  ├─ Identity threshold: >=90%
   │  ├─ Coverage threshold: >=60%
   │  └─ Plasmid replicons detected:
   │     ├─ rep7a    -- Most prevalent; South Korea, Brazil, UK, Jersey, Vietnam, Australia
   │     ├─ rep7b    -- Restricted to Vietnamese clinical strains (SDH.B1, SDH.B3)
   │     ├─ rep19c   -- Associated with lnu(A) carriage
   │     └─ rep21    -- Broad distribution; associated with lnu(A)
   │
   └─ sal(E) gene analysis
      ├─ ROTIFER v1.0
      │  ├─ Extract sal(E) sequences from all 36 genomes
      │  ├─ Extract ~2 kb flanking regions on each side
      │  └─> 36 sal(E)-containing sequences with flanking context
      │
      ├─ Pfam database annotation (v35.0)
      │  ├─ Annotate flanking ORFs with conserved domains
      │  ├─ Upstream genes:   iscS/nifS, mnmA, TPR proteins, recD2
      │  ├─ Downstream genes: alaS, ruvX, aspS, ssrS (6S RNA)
      │  └─> Conserved synteny across all 36 genomes -- single ancestral acquisition
      │
      ├─ MAFFT v7.525 (multiple sequence alignment)
      │  ├─ 36 sal(E) protein sequences aligned
      │  ├─ 543 amino acid positions
      │  ├─ Mean amino acid identity: >92%
      │  ├─ 531 constant sites (97.8%)
      │  └─> 7 amino acid substitutions + 2 deletions across dataset
      │     at positions: Cys179, His229, Arg312, His386,
      │     His431, Asp450, Ile473 + deletions at 60 and 505
      │
      ├─ AlphaFold3 (web server: https://alphafoldserver.com)
      │  ├─ In silico structural prediction of Sal(E)
      │  ├─ pTM = 0.8 (high confidence)
      │  ├─ Structurally identical to Sal(B) cryo-EM structure
      │  └─> Supports canonical ABC-F fold and ribosomal target-protection mechanism
      │
      └─ IQ-TREE v3.1.2 (gene-specific phylogeny)
         ├─ Input: sal_E_aligned.faa (36 sequences, 543 aa)
         ├─ Model: Q.YEAST (selected by ModelFinder)
         ├─ Bootstrap: 1,000 UFBoot
         ├─ Constant sites: 531/543 (97.8%)
         ├─ Parsimony-informative sites: 3
         └─> Highly conserved gene tree consistent with
             single ancestral acquisition and strict vertical inheritance

6. Core Genome Analysis
   ├─ Input: core_gene_alignment_filtered.aln from Panaroo
   │         (2,348 core genes; 1,677,386 total positions)
   │
   ├─ Gubbins v3.4.2
   │  ├─ Recombination detection and masking
   │  ├─ Outgroup: GCF_002442935.1 (J11 reference strain)
   │  ├─ Tree builder: FastTree
   │  ├─ Total SNPs:              32,454
   │  ├─ SNPs in recombination:    1,908  (5.88%)
   │  ├─ Recombination blocks:       211
   │  ├─ Mean r/m per branch:       0.100
   │  ├─ Mean rho/theta:            0.015
   │  └─> nepalensis.filtered_polymorphic_sites.fasta → input for IQ-TREE
   │
   └─ ClonalFrameML v1.20
      ├─ Input: Gubbins final tree + core alignment
      ├─ emsim: 100 bootstrap simulations
      ├─ R/theta:   3.97 × 10⁻⁵
      ├─ delta:     97.58 bp  (mean recombinant tract length)
      ├─ nu:        2.09 × 10⁻⁴  (mean import divergence)
      ├─ r/m:       8.09 × 10⁻⁷  (= R/theta × delta × nu)
      └─> Confirms predominantly clonal population structure

7. Phylogenetic Inference
   ├─ SNP-sites v2.5.1
   │  ├─ Input: Gubbins filtered alignment (1,677,386 positions)
   │  ├─ Extracts variable positions only
   │  └─> 20,043 SNP positions retained
   │
   ├─ SNP-dists v0.9.0
   │  ├─ Input: Gubbins filtered alignment
   │  └─> 36 × 36 pairwise SNP distance matrix
   │     ├─ COLB vs AM1:                    1 cgSNP (clonal)
   │     ├─ COLB/AM1 vs Korean strains:  ≥1,648 cgSNPs
   │     └─ COLB/AM1 vs all others:      ≥3,491 cgSNPs
   │
   └─ IQ-TREE v3.1.2
      ├─ Input: nepalensis.filtered_polymorphic_sites.fasta
      │         (1,677,386 positions; 20,043 variable; 8,739 parsimony-informative)
      ├─ Model: GTR+G (fixed)
      ├─ Constant site correction: -fconst 590159,246851,316673,503660
      ├─ Bootstrap: 1,000 ultrafast (UFBoot2)
      ├─ Best log-likelihood: -2,471,604.131
      └─> nepalensis.treefile → visualized in iTOL v6

8. Output & Results
   ├─ Phylogenetic trees
   │  ├─ nepalensis.treefile         (IQ-TREE ML tree -- Figure 1B)
   │  ├─ nepalensis.final_tree.tre   (Gubbins recombination-corrected -- Figure 2A)
   │  └─ sal_E.treefile              (sal(E) gene tree -- Figure 3)
   │
   ├─ Distance matrices
   │  ├─ nepalensis.snp_distance_matrix.tsv  (cgSNP pairwise)
   │  └─ fastani_results.txt.matrix          (ANI all-vs-all)
   │
   ├─ AMR profiles
   │  ├─ nepalensis_complete_amr_table_v2.csv  (complete metadata + AMR)
   │  └─ recombination_summary.csv             (Gubbins + ClonalFrameML)
   │
   └─ sal(E) genetic context
      ├─ Conserved synteny documented across all 36 genomes
      ├─ Upstream:   iscS/nifS -> mnmA -> TPR -> recD2
      ├─ Downstream: alaS -> ruvX -> aspS -> ssrS
      └─ AlphaFold3 structural model (pTM = 0.8)
         consistent with ABC-F ribosomal target-protection mechanism

```

---

## Requirements

### Software

All tools are managed via conda/mamba environments:

```bash
# Genome download
mamba create -n download        -c bioconda ncbi-datasets-cli

# Quality assessment
mamba create -n checkm2         -c bioconda checkm2
mamba create -n quast           -c bioconda quast=5.2.0

# Genome annotation
mamba create -n bakta           -c bioconda bakta=1.12.0

# ANI analysis
mamba create -n fastani         -c bioconda fastani

# Pan-genome analysis
mamba create -n panaroo         -c bioconda panaroo

# Recombination detection and parameters
mamba create -n gubbins         -c bioconda gubbins
mamba create -n clonalframe     -c bioconda clonalframeml

# SNP extraction and distances
mamba create -n snpsites        -c bioconda snp-sites snp-dists

# Phylogenetic inference
mamba create -n iqtree          -c bioconda iqtree

# AMR analysis (--plus flag enables heavy metal / stress gene detection)
mamba create -n amrfinder       -c bioconda ncbi-amrfinderplus
mamba create -n card_new        -c bioconda rgi

# Plasmid replicon typing
mamba create -n plasmidfinder   -c bioconda plasmidfinder

# Virulence detection
mamba create -n abricate        -c bioconda abricate
mamba create -n virulencefinder -c bioconda virulencefinder

# sal(E) genetic context and structural analysis
mamba create -n alignment       -c bioconda mafft iqtree
mamba create -n rotifer         -c bioconda rotifer

# Note: AlphaFold3 accessed via web server https://alphafoldserver.com
# No local installation required
```

### Databases

| Database           | Version           | Path                                              |
| ------------------ | ----------------- | ------------------------------------------------- |
| Bakta full DB      | v6.0 (2025-02-24) | `/databases/bakta_db/db`                          |
| CheckM2            | 2026-05-15.1      | `/databases/CheckM2_database/uniref100.KO.1.dmnd` |
| AMRFinderPlus      | 2026-05-15.1      | auto-managed                                      |
| CARD               | 3.2.7             | `/databases/card/card.json`                       |
| PlasmidFinder DB   | 2024-01-24        | `/databases/plasmidfinder_db`                     |
| VirulenceFinder DB | v2.0.1            | `/databases/virulencefinder_db`                   |
| Pfam               | v35.0             | auto-managed / `/databases/pfam/Pfam-A.hmm`       |

---

## Step-by-Step Pipeline

### Step 1 — Download Genomes

```bash
conda activate download
datasets download genome taxon "Staphylococcus nepalensis"
unzip ncbi_dataset.zip -d nepalensis_ncbi_genomes
# Result: 88 genomes
```

### Step 2 — Deduplication

```bash
# Run from: nepalensis_ncbi_genomes/ncbi_dataset/data/
# Groups by BioSample, keeps highest version, prefers GCF over GCA
# Result: 46 unique genomes, 42 duplicates moved to duplicates/

bash ../../scripts/dedup.sh
# or manually:

jq -r '.accession + "\t" + (.assemblyInfo.biosample.accession // "NO_BIOSAMPLE")' \
    assembly_data_report.jsonl > /tmp/assembly_biosample.txt

# Implementation: see scripts/dedup.sh for full logic
# - Groups by BioSample
# - Prefers GCF over GCA
# - Keeps highest version
# - Moves duplicates to duplicates/ folder
```

### Step 3 — Export Unique FASTAs

```bash
mkdir -p exported_fastas
jq -r '.accession' assembly_data_report.jsonl | while read acc; do
    fna=$(find "$acc" -name "*.fna" | head -1)
    [ -f "$fna" ] && cp "$fna" "exported_fastas/${acc}.fna"
done
# Result: 46 .fna files
```

### Step 4 — Create Metadata Table

```bash
# Extracts: accession, biosample, taxid, organism, strain,
# year, country, isolation_source, assembly_level from NCBI JSON
jq -r '... | @csv' assembly_data_report.jsonl > genome_tracking.csv
```

### Step 5 — Quality Assessment

```bash
# CheckM2
conda activate checkm2
checkm2 database --setdblocation /databases/CheckM2_database/uniref100.KO.1.dmnd
checkm2 predict \
    --input exported_fastas/ \
    --output-directory checkm2_results \
    --extension fna \
    --threads 10

# QUAST
conda activate quast
quast exported_fastas/*.fna \
    --output-dir quast_results \
    --threads 10
```

### Step 6 — Quality Filtering (Python)

```python
# filter_genomes.py
# Thresholds:
#   Completeness >= 95%  (CheckM2)
#   Contamination <= 5%  (CheckM2)
#   Total contigs <= 100 (QUAST)
#   Total N's <= 500     (QUAST)
#   Exclude suppressed   (NCBI)
# Result: 36 high-quality genomes

import pandas as pd
import numpy as np

df = pd.read_csv("nepalensis_complete_table.csv")
df["Total_Ns"] = (df["# N's per 100 kbp"] / 100000) * df["Total length"]
df["Total_Ns"] = df["Total_Ns"].round(0).astype(int)

mask = (
    (df["Completeness"] >= 95) &
    (df["Contamination"] <= 5) &
    (df["Total_Contigs"] <= 100) &
    (df["Total_Ns"] <= 500) &
    (df["version_status"] != "suppressed")
)
df[mask].to_csv("nepalensis_highquality.csv", index=False)
df[~mask].to_csv("nepalensis_lowquality.csv", index=False)
```

### Step 7 — Genome Annotation (Bakta)

```bash
conda activate bakta
bakta_db --setdblocation /databases/bakta_db/db

annotate_genome() {
    fna=$1
    accession=$(basename "$fna" .fna)
    bakta \
        --db /databases/bakta_db/db \
        --output bakta_results/${accession} \
        --prefix ${accession} \
        --threads 4 \
        --force \
        --skip-sorf \
        "$fna"
}
export -f annotate_genome
parallel -j 3 annotate_genome ::: hq_fastas/*.fna

# Note: --skip-sorf required due to DIAMOND segfault with bakta 1.12.0 + Python 3.13
# Workaround: Use bakta 1.12.0 with Python 3.12 or upgrade bakta to v1.13.0+
# See: https://github.com/oschwengers/bakta/issues/XXX
# Result: 36 annotated genomes (.gff3, .gbff, .faa, .tsv, .png)
```

### Step 8 — Average Nucleotide Identity (FastANI)

```bash
conda activate fastani

ls hq_fastas/*.fna > fastani_genome_list.txt

# All vs all
fastANI \
    --ql fastani_genome_list.txt \
    --rl fastani_genome_list.txt \
    --output fastani_results.txt \
    --threads 10 \
    --matrix

# All vs reference genome (GCF_002442935.1 = NCBI reference strain J11)
fastANI \
    --ql fastani_genome_list.txt \
    --ref hq_fastas/GCF_002442935.1.fna \
    --output fastani_vs_reference.txt \
    --threads 10
```

### Step 9 — Pangenome Analysis (Panaroo)

```bash
conda activate panaroo

panaroo \
    -i bakta_results/*/*.gff3 \
    -o panaroo_results \
    --clean-mode strict \
    -a core \
    --core_threshold 0.98 \
    -t 20

# Results:
#   Core genes (99-100%):     2,348
#   Soft core genes (95-99%):    40
#   Shell genes (15-95%):       485
#   Cloud genes (<15%):       1,880
#   Total pangenome:          4,753
#   Pangenome type:           Open
```

### Step 10 — Recombination Detection (Gubbins)

```bash
conda activate gubbins

run_gubbins.py \
    --prefix nepalensis \
    --tree-builder fasttree \
    --threads 20 \
    --outgroup GCF_002442935.1 \
    panaroo_results/core_gene_alignment_filtered.aln

# Results:
#   Total SNPs:                32,454
#   SNPs in recombination:      1,908 (5.88%)
#   Recombination blocks:         211
#   Mean r/m:                   0.100
#   Mean rho/theta:             0.015
#
# Key output: nepalensis.filtered_polymorphic_sites.fasta
# Used as input for SNP-sites and IQ-TREE
```

### Step 11 — Recombination Parameters (ClonalFrameML)

```bash
conda activate clonalframe

ClonalFrameML \
    gubbins_results/nepalensis.final_tree.tre \
    panaroo_results/core_gene_alignment_filtered.aln \
    clonalframe_results/nepalensis \
    -emsim 100 \
    -num_threads 20 \
    -show_progress true

# Results:
#   R/theta:  3.97e-05
#   delta:    97.58 bp  (mean recombinant tract length)
#   nu:       2.09e-04  (mean divergence of imports)
#   r/m:      8.09e-07
#
# Formula: r/m = R/theta x delta x nu
```

### Step 12 — SNP Extraction (SNP-sites)

```bash
conda activate snpsites

# Extract SNP-only alignment from Gubbins filtered alignment
snp-sites \
    -o gubbins_results/nepalensis.snps.aln \
    gubbins_results/nepalensis.filtered_polymorphic_sites.fasta

# Count constant sites for IQ-TREE branch length correction
snp-sites -C panaroo_results/core_gene_alignment_filtered.aln
# Output: 590159,246851,316673,503660 (A,C,G,T)
# These values are passed to IQ-TREE with -fconst
```

### Step 13 — Phylogenetic Analysis (IQ-TREE)

```bash
conda activate iqtree

# Input: Gubbins recombination-filtered alignment (1,677,386 positions)
# Constant sites correction applied with -fconst
iqtree \
    -s gubbins_results/nepalensis.filtered_polymorphic_sites.fasta \
    -m GTR+G \
    -fconst 590159,246851,316673,503660 \
    -bb 1000 \
    -nt 20 \
    -pre iqtree_results/nepalensis

# Results:
#   Input:               Gubbins recombination-filtered alignment
#   Total positions:     1,677,386
#   Variable sites:         20,043
#   Parsimony-informative:   8,739
#   Constant sites:      1,657,343 (corrected via -fconst)
#   Bootstrap:           1,000 UFBoot2
#   Model:               GTR+G
#   Best log-likelihood: -2,471,604.131
```

### Step 14 — Pairwise SNP Distance Matrix (SNP-dists)

```bash
conda activate snpsites

snp-dists \
    gubbins_results/nepalensis.filtered_polymorphic_sites.fasta \
    > gubbins_results/nepalensis.snp_distance_matrix.tsv

# Result: 36x36 pairwise SNP distance matrix
# Key distances:
#   COLB vs AM1:                   1 cgSNP (clonal)
#   COLB/AM1 vs Korean strains: >=1,648 cgSNPs
#   COLB/AM1 vs all others:     >=3,491 cgSNPs
# Identical genome pairs detected:
#   GCF_900458695.1 = GCF_014635045.1
#   GCF_051590305.1 = GCF_051590285.1
```

### Step 15 — Antimicrobial Resistance Analysis

#### AMRFinderPlus

```bash
conda activate amrfinder

mkdir -p amrfinder_results
for faa in bakta_results/*/*.faa; do
    accession=$(basename "$faa" .faa)
    [[ "$accession" == *.hypotheticals ]] && continue
    amrfinder \
        -p "$faa" \
        --organism Staphylococcus_aureus \
        --plus \
        --output amrfinder_results/${accession}.tsv \
        --threads 4
done

# Note: --plus flag required to detect heavy metal and stress resistance genes

# Core resistome (100% genomes):
#   sal(E)  -- Lincosamide/Pleuromutilin/Streptogramin A
#
# Accessory AMR genes:
#   fosB/fosB4         21/36 (58.3%) -- Fosfomycin
#   tet(K)              9/36 (25.0%) -- Tetracycline
#   mph(C)              6/36 (16.7%) -- Macrolide
#   qacG                5/36 (13.9%) -- Quaternary ammonium (biocide)
#   str                 4/36 (11.1%) -- Aminoglycoside
#   lnu(A)/lnu(A)'      3/36  (8.3%) -- Lincosamide
#   catA                3/36  (8.3%) -- Phenicol
#   aac(6')-Ie/aph(2'') 2/36  (5.6%) -- Aminoglycoside
#   mecA/mecI/mecR1     2/36  (5.6%) -- Beta-lactam/MRSA
#   dfrE                2/36  (5.6%) -- Trimethoprim
#   erm(B)              2/36  (5.6%) -- Macrolide
#   dfrG                1/36  (2.8%) -- Trimethoprim
#   blaI/blaPC1/blaR1   1/36  (2.8%) -- Beta-lactam regulators
#
# Heavy metal & stress resistance genes (detected via --plus flag):
#   arsR/arsB/arsC  -- Arsenic resistance operon (widely distributed)
#   cadD            -- Cadmium tolerance
#   mco             -- Copper tolerance / oxidative stress
#   merA/merB/merT  -- Mercury resistance operon
```

#### CARD (RGI)

```bash
conda activate card_new
rgi load -i /databases/card/card.json --local

mkdir -p card_results
for faa in bakta_results/*/*.faa; do
    accession=$(basename "$faa" .faa)
    [[ "$accession" == *.hypotheticals ]] && continue
    rgi main \
        -i "$faa" \
        -o card_results/${accession} \
        -t protein \
        -a BLAST \
        --local \
        --clean \
        -n 4
done
# Reporting: Perfect and Strict hits only
# Confirmed: sal(E), tet(K), mecA, lnu(A), mph(C), qacG
```

#### Plasmid Replicon Typing (PlasmidFinder)

```bash
conda activate plasmidfinder

# Download database if not available
cd /databases/plasmidfinder_db
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git .

mkdir -p plasmidfinder_results
for fna in hq_fastas/*.fna; do
    accession=$(basename "$fna" .fna)
    mkdir -p plasmidfinder_results/${accession}
    plasmidfinder.py \
        -i "$fna" \
        -o plasmidfinder_results/${accession} \
        -p /databases/plasmidfinder_db \
        -l 0.60 \
        -t 0.90
done

# Plasmid replicons detected:
#   rep7a    -- Most prevalent; South Korea, Brazil, UK, Jersey, Vietnam, Australia
#   rep7b    -- Restricted to Vietnamese clinical strains (SDH.B1, SDH.B3)
#   rep19c   -- Associated with lnu(A) carriage
#   rep21    -- Broad distribution; associated with lnu(A)
```

#### sal(E) Genetic Context Analysis

##### Extract sal(E) Flanking Sequences (ROTIFER)

```bash
conda activate rotifer

mkdir -p gene_extraction/sale_context

rotifer extract \
    --input hq_fastas/ \
    --gene "sal(E)" \
    --flank 2000 \
    --output gene_extraction/sale_context/

# Result: 36 sal(E)-containing sequences with ~2 kb flanking regions
# on each side for genetic context analysis
```

##### Annotate Flanking Sequences with Pfam

```bash
conda activate alignment

mkdir -p gene_extraction/sale_context/pfam_annotations

hmmscan \
    --domtblout gene_extraction/sale_context/pfam_annotations/flanking.domtbl \
    /databases/pfam/Pfam-A.hmm \
    gene_extraction/sale_context/*_flanking.faa

python3 scripts/extract_orfs_flanking.py \
    gene_extraction/sale_context/ \
    gene_extraction/sale_context/pfam_annotations/

# Results:
#   Upstream genes:   iscS/nifS, mnmA, TPR proteins, recD2
#   Downstream genes: alaS, ruvX, aspS, ssrS (6S RNA)
#   Conserved synteny across all 36 genomes
#   Evidence of single ancestral acquisition
```

### Step 16 — Virulence Factor Detection

#### VirulenceFinder

```bash
conda activate virulencefinder

mkdir -p /databases/virulencefinder_db
cd /databases/virulencefinder_db
git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git .
kma index -i s.aureus_toxin.fsa -o s.aureus_toxin
kma index -i s.aureus_exoenzyme.fsa -o s.aureus_exoenzyme
kma index -i s.aureus_hostimm.fsa -o s.aureus_hostimm

mkdir -p virulencefinder_results
for fna in hq_fastas/*.fna; do
    accession=$(basename "$fna" .fna)
    mkdir -p virulencefinder_results/${accession}
    virulencefinder.py \
        -i "$fna" \
        -o virulencefinder_results/${accession} \
        -p /databases/virulencefinder_db \
        -d s.aureus_toxin,s.aureus_exoenzyme,s.aureus_hostimm \
        -l 0.6 \
        -t 0.9 \
        -x
done
# Result: No S. aureus virulence genes detected in any genome
```

#### Abricate (VFDB)

```bash
conda activate abricate

abricate \
    --db vfdb \
    --minid 80 \
    --mincov 60 \
    hq_fastas/*.fna \
    > abricate_results/vfdb_results.tsv

# Result: No virulence genes detected
# Interpretation: S. nepalensis lacks known S. aureus virulence factors,
# consistent with its non-pathogenic CoNS nature
```

### Step 17 — Gene-specific Phylogenetic Analysis (sal(E))

```bash
# Extract sal(E) sequences from all 36 genomes
conda activate base
python3 scripts/extract_sale.py
# Result: 36 sal(E) sequences extracted

# Multiple sequence alignment
conda activate alignment
mafft --auto sal_E_all.faa > sal_E_aligned.faa

# Phylogenetic tree
iqtree \
    -s sal_E_aligned.faa \
    -m TEST \
    -bb 1000 \
    -T AUTO \
    -pre sal_E_tree/sal_E

# Results:
#   Best model:              Q.YEAST
#   Total sites:             543 aa
#   Constant sites:          531 (97.8%)
#   Parsimony-informative:   3
#   Mean amino acid identity: >92%
#   Conclusion: sal(E) is highly conserved -- single ancestral acquisition
#               consistent with strict vertical inheritance
```

### Step 18 — Merge All Results (Python)

Provided scripts in `scripts/` directory:

```python
# Scripts available:
# merge_fastani.py         -- adds FastANI vs reference to high-quality table
# merge_amr_tables_v2.py   -- adds AMRFinder + CARD + PlasmidFinder results
# recombination_summary.py -- summarizes Gubbins + ClonalFrameML stats
# copy_hq_fastas.py        -- copies 36 HQ FASTAs to separate folder
# merge_sale_context.py    -- adds sal(E) genetic context annotations
```

### Step 19 — Transfer Results to Local Machine

```bash
# Run from local machine
mkdir -p /home/johana/Documentos/nepalensis

scp -r zbecerra@davinci.icb.usp.br:/path/to/export_for_johana \
    /home/johana/Documentos/nepalensis/
scp -r zbecerra@davinci.icb.usp.br:/path/to/gubbins_results \
    /home/johana/Documentos/nepalensis/
scp -r zbecerra@davinci.icb.usp.br:/path/to/iqtree_results \
    /home/johana/Documentos/nepalensis/
scp -r zbecerra@davinci.icb.usp.br:/path/to/clonalframe_results \
    /home/johana/Documentos/nepalensis/
scp -r zbecerra@davinci.icb.usp.br:/path/to/bakta_results \
    /home/johana/Documentos/nepalensis/
scp -r zbecerra@davinci.icb.usp.br:/path/to/amrfinder_results \
    /home/johana/Documentos/nepalensis/
scp -r zbecerra@davinci.icb.usp.br:/path/to/card_results \
    /home/johana/Documentos/nepalensis/
scp -r zbecerra@davinci.icb.usp.br:/path/to/plasmidfinder_results \
    /home/johana/Documentos/nepalensis/
scp -r zbecerra@davinci.icb.usp.br:/path/to/gene_extraction \
    /home/johana/Documentos/nepalensis/
```

---

## Output Files

```
projects/nepalensis/
├── export_for_johana/
│   ├── exported_fastas/                     # 46 unique FASTAs
│   ├── hq_fastas/                           # 36 high-quality FASTAs
│   ├── genome_tracking.csv                  # NCBI metadata (46 genomes)
│   ├── nepalensis_complete_table.csv        # metadata + CheckM2 + QUAST
│   ├── nepalensis_highquality.csv           # filtered table (36 genomes)
│   ├── nepalensis_lowquality.csv            # excluded genomes + reasons
│   ├── nepalensis_highquality_fastani.csv   # table + ANI vs reference
│   ├── nepalensis_complete_amr_table_v2.csv # complete table + AMR
│   ├── recombination_summary.csv            # Gubbins + ClonalFrameML
│   └── recombination_per_branch.csv         # per-branch statistics
│
├── checkm2_results/
│   └── quality_report.tsv
│
├── quast_results/
│   ├── transposed_report.tsv
│   └── report.html
│
├── bakta_results/
│   └── <accession>/
│       ├── <accession>.gff3
│       ├── <accession>.gbff
│       ├── <accession>.faa
│       ├── <accession>.ffn
│       ├── <accession>.tsv
│       └── <accession>.png
│
├── panaroo_results/
│   ├── gene_presence_absence.csv
│   ├── core_gene_alignment.aln
│   ├── core_gene_alignment_filtered.aln
│   └── summary_statistics.txt
│
├── gubbins_results/
│   ├── nepalensis.final_tree.tre
│   ├── nepalensis.recombination_predictions.gff
│   ├── nepalensis.per_branch_statistics.csv
│   ├── nepalensis.filtered_polymorphic_sites.fasta
│   ├── nepalensis.snps.aln
│   └── nepalensis.snp_distance_matrix.tsv
│
├── clonalframe_results/
│   ├── nepalensis.em.txt
│   ├── nepalensis.emsim.txt
│   ├── nepalensis.labelled_tree.newick
│   ├── nepalensis.importation_status.txt
│   └── nepalensis.ML_sequence.fasta
│
├── iqtree_results/
│   ├── nepalensis.treefile
│   ├── nepalensis.contree
│   ├── nepalensis.iqtree
│   └── nepalensis.mldist
│
├── amrfinder_results/
│   ├── <accession>.tsv
│   └── amrfinder_all_genomes.tsv
│
├── card_results/
│   ├── <accession>.txt
│   └── card_all_genomes.tsv
│
├── plasmidfinder_results/              
│   └── <accession>/
│       └── results_tab.tsv
│
├── virulencefinder_results/
│   └── <accession>/
│       └── data.json
│
├── abricate_results/
│   └── vfdb_results.tsv
│
├── gene_extraction/
│   ├── sal_E_all.faa
│   ├── sal_E_aligned.faa
│   ├── sale_context/
│   │   ├── <accession>_sale.fna
│   │   ├── <accession>_flanking.fna
│   │   ├── <accession>_flanking.faa
│   │   └── pfam_annotations/
│   │       ├── flanking.domtbl
│   │       └── genetic_context_summary.csv
│   └── sal_E_tree/
│       ├── sal_E.treefile
│       └── sal_E.iqtree
│
├── fastani_results.txt
├── fastani_results.txt.matrix
├── fastani_vs_reference.txt
└── fastani_genome_list.txt
```

---

## Key Results Summary

| Analysis                           | Result                                   |
| ---------------------------------- | ---------------------------------------- |
| Total genomes downloaded           | 88                                       |
| Unique genomes after deduplication | 46                                       |
| High-quality genomes               | 36                                       |
| Core genome size                   | 2,348 genes                              |
| Soft core genes                    | 40 genes                                 |
| Shell genes                        | 485 genes                                |
| Cloud genes                        | 1,880 genes                              |
| Total pangenome                    | 4,753 genes                              |
| Pangenome type                     | Open                                     |
| Gubbins SNPs (total)               | 32,454                                   |
| SNPs used for phylogeny (IQ-TREE)  | 20,043                                   |
| Recombination blocks (Gubbins)     | 211                                      |
| SNPs from recombination            | 1,908 (5.88%)                            |
| Gubbins r/m                        | 0.100                                    |
| ClonalFrameML R/theta              | 3.97 x 10⁻⁵                              |
| ClonalFrameML delta                | 97.58 bp                                 |
| ClonalFrameML nu                   | 2.09 x 10⁻⁴                              |
| ClonalFrameML r/m                  | 8.09 x 10⁻⁷                              |
| Core resistome (100% genomes)      | sal(E), arr, bla                         |
| Most prevalent accessory AMR       | fosB/fosB4 (58.3%), tet(K) (25.0%)       |
| Heavy metal resistance             | arsR/B/C, cadD, mco, merA/B/T           |
| Biocide resistance                 | qacG (13.9%) -- clinical isolates only   |
| Most prevalent plasmid replicon    | rep7a                                    |
| MDR-associated replicon            | rep7b -- Vietnamese clinical strains only|
| sal(E) conservation                | 97.8% (531/543 amino acid sites)         |
| sal(E) genetic context             | Highly conserved across all genomes      |
| Virulence genes detected           | None                                     |

---

## sal(E) Genetic Context Analysis

The sal(E) gene and its flanking sequences were extracted from all 36 genomes using
ROTIFER (Rapid Open-source Tools and Infrastructure for Data Exploration and Research;
https://github.com/leepbioinfo/rotifer). Approximately 2 kb regions flanking sal(E) on
both sides were retrieved to determine the genetic organization and evolutionary
conservation of the locus.

The extracted flanking sequences were annotated using the Pfam database to identify
conserved protein domains in adjacent coding regions. This analysis revealed:

- **Genetic organization**: Consistent synteny across all S. nepalensis isolates
- **Upstream genes**: iscS/nifS → mnmA → TPR proteins → recD2
- **Downstream genes**: alaS → ruvX → aspS → ssrS (6S RNA)
- **Evolutionary conservation**: Evidence of single ancestral acquisition with minimal
  subsequent modification
- **Functional context**: Association with tRNA modification, translation fidelity,
  and DNA repair pathways
- **Structural prediction**: AlphaFold3 model (pTM = 0.8) confirms ABC-F fold
  consistent with ribosomal target-protection mechanism

---

## Recombination Parameters

```
r/m = R/theta × delta × nu
r/m = 3.97e-05 × 97.58 × 2.09e-04
r/m = 8.09e-07

Interpretation: S. nepalensis is a predominantly clonal species with
very low recombination rates. Short recombinant tracts (97.58 bp) and
low import divergence (2.09e-04) suggest within-species recombination
rather than importation from distantly related organisms, consistent
with other coagulase-negative staphylococci.
```

---

## Phylogenetic Analysis Notes

The maximum likelihood phylogeny was inferred using IQ-TREE v3.1.2 with the GTR+G
substitution model. The input alignment consisted of the full recombination-filtered
core-genome alignment produced by Gubbins (1,677,386 positions; 20,043 variable sites;
8,739 parsimony-informative sites). Constant sites (A=590,159; C=246,851; G=316,673;
T=503,660) were corrected using the -fconst option. Bootstrap support was assessed with
1,000 ultrafast bootstrap replicates (UFBoot2).

---

## Virulence Factor Analysis

No virulence genes were detected in any of the 36 S. nepalensis genomes using:

- **VirulenceFinder** (databases: s.aureus_toxin, s.aureus_exoenzyme, s.aureus_hostimm;
  identity ≥ 90%, coverage ≥ 60%)
- **Abricate** with VFDB (identity ≥ 80%, coverage ≥ 60%)

This confirms the non-pathogenic nature of S. nepalensis and its distinction from
pathogenic staphylococci such as *S. aureus*, consistent with its classification as a
coagulase-negative commensal species.

---

## Software Versions

| Tool | Version | Reference |
|---|---|---|
| NCBI Datasets CLI | 18.29.1 | NCBI |
| CheckM2 | 1.1.0 | Chklovski et al. 2023 |
| QUAST | 5.2.0 | Gurevich et al. 2013 |
| Bakta | 1.12.0 | Schwengers et al. 2021 |
| FastANI | 1.34 | Jain et al. 2018 |
| Panaroo | 1.3.4 | Tonkin-Hill et al. 2020 |
| Gubbins | 3.4.2 | Croucher et al. 2015 |
| ClonalFrameML | 1.20 | Didelot & Wilson 2015 |
| SNP-sites | 2.5.1 | Page et al. 2016 |
| SNP-dists | 0.9.0 | Seemann 2018 |
| IQ-TREE | 3.1.2 | Minh et al. 2020 |
| AMRFinderPlus | 4.2.7 | Feldgarden et al. 2021 |
| RGI/CARD | - (DB: 3.2.7) | Alcock et al. 2023 |
| PlasmidFinder | 2.2.0 | Carattoli et al. 2014 |
| ROTIFER | 1.0 | https://github.com/leepbioinfo/rotifer |
| Pfam | v35.0 | Mistry et al. 2021 |
| VirulenceFinder | 2.0.1 | Joensen et al. 2014 |
| Abricate | latest | Seemann 2020 |
| MAFFT | 7.525 | Katoh & Standley 2013 |
| AlphaFold3 | web server | Abramson et al. 2024 |

---

## References

- Abramson et al. (2024) Accurate structure prediction of biomolecular interactions with AlphaFold3. Nature
- Alcock et al. (2023) CARD 2023: Expanded antibiotic resistance gene database. Nucleic Acids Research
- Carattoli et al. (2014) PlasmidFinder and pMLST: in silico detection and typing of plasmids. Antimicrobial Agents and Chemotherapy
- Chklovski et al. (2023) CheckM2: A rapid, scalable and accurate assessment of microbial genome quality using machine learning. Nature Methods
- Croucher et al. (2015) Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. Nucleic Acids Research
- Didelot & Wilson (2015) ClonalFrameML: efficient inference of recombination in whole bacterial genomes. PLoS Computational Biology
- Feldgarden et al. (2021) AMRFinderPlus: a database and bioinformatics pipeline for antibiotic resistance determinants. Journal of Antimicrobial Chemotherapy
- Gurevich et al. (2013) QUAST: quality assessment tool for genome assemblies. Bioinformatics
- Jain et al. (2018) High-throughput ANI analysis of 90K prokaryotic genomes reveals sharp species boundaries. Nature Communications
- Joensen et al. (2014) Real-time whole-genome sequencing for routine typing, surveillance, and outbreak detection of verotoxigenic Escherichia coli. Journal of Clinical Microbiology
- Katoh & Standley (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution
- Minh et al. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution
- Mistry et al. (2021) Pfam: the protein families database in 2021. Nucleic Acids Research
- Page et al. (2016) SNP-sites: rapid efficient extraction of SNPs from multi-FASTA alignments. Microbial Genomics
- Schwengers et al. (2021) Bakta: rapid and standardized annotation of bacterial genomes via long-read sequencing. Microbial Genomics
- Seemann (2020) Abricate: mass screening of contigs for antimicrobial and virulence genes. https://github.com/tseemann/abricate
- Tonkin-Hill et al. (2020) Producing polished prokaryotic pangenomes with the Panaroo pipeline. Genome Biology

---

## Citation

If you use this pipeline, please cite:

> Becerra et al. (2026). Phylogeographic analysis of Staphylococcus nepalensis reveals global occurrence of antimicrobial-resistant lineages carrying the intrinsic sal(E) resistance gene. [Journal name and volume to be added]

---

## Troubleshooting

### Common Issues

**Bakta DIAMOND segfault (Python 3.13 + Bakta 1.12.0)**
- **Solution 1**: Use `--skip-sorf` flag (already included in Step 7)
- **Solution 2**: Downgrade to Python 3.12 in the bakta environment
- **Solution 3**: Upgrade to Bakta v1.13.0+ (recommended for new installations)

**Database path errors**
- Ensure conda environments have correct `--setdblocation` paths
- Verify database files exist: `ls /databases/bakta_db/db`, `ls /databases/CheckM2_database/`
- For custom paths, update paths in all scripts before running

**Memory issues with Panaroo or Gubbins**
- Reduce thread count (`-t`) parameter
- Process genomes in smaller batches if necessary
- Ensure adequate disk space for temporary files (>50GB recommended)

### Contact & Support

For pipeline-specific issues, see `scripts/` directory for helper scripts.
For tool-specific help, consult official documentation:
- Bakta: https://github.com/oschwengers/bakta
- Panaroo: https://github.com/gtonkinhill/panaroo
- Gubbins: https://github.com/sanger-pathogens/gubbins
- IQ-TREE: http://www.iqtree.org/

---

## License

MIT License
