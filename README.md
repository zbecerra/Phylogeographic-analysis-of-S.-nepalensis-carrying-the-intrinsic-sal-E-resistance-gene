# Phylogeographic analysis of Staphylococcus nepalensis reveals global occurrence of antimicrobial-resistant lineages carrying the intrinsic sal(E) resistance gene*

## Overview
This repository contains the complete bioinformatics pipeline used for the comparative genomic analysis of *Staphylococcus nepalensis*, including genome quality assessment, annotation, pangenome analysis, recombination detection, phylogenetics, antimicrobial resistance profiling, virulence factor detection, and gene-specific analyses.

---

## Pipeline Overview

```
1. NCBI Genome Download
   └─> 88 genomes

2. Deduplication by BioSample
   └─> 46 unique genomes (42 duplicates removed)

3. Quality Filtering (CheckM2 + QUAST)
   ├─ Completeness >= 95%
   ├─ Contamination <= 5%
   ├─ Total contigs <= 100
   └─ Total N's <= 500
   └─> 36 high-quality genomes

4. Parallel Processing
   ├─ FastANI (ANI matrix)
   ├─ Bakta (Annotation)
   ├─ Panaroo (Pangenome analysis)
   └─ VirulenceFinder + Abricate VFDB (Virulence detection)

5. AMR Analysis
   ├─ AMRFinderPlus + CARD
   ├─ sal(E) gene analysis (Core resistome: 100% genomes)
   │  ├─ ROTIFER (Rapid Open-source Tools and Infrastructure
   │  │           for Data Exploration and Research)
   │  │  └─ Extract nucleotide sequences flanking sal(E)
   │  │  └─ Determine genetic context and organization
   │  │
   │  └─ Pfam database annotation
   │     └─ Annotate flanking sequences and coding regions
   │
   └─> Accessory resistance genes (fosB 64%, tet(K) 25%, 
                                   mph(C) 17%, mecA 6%)

6. Core Genome Analysis
   ├─ Core gene alignment from Panaroo
   ├─ Gubbins (Recombination detection)
   ├─ ClonalFrameML (Recombination parameters)
   └─> 32,454 total SNPs

7. Phylogenetic Inference
   ├─ SNP-sites (Extract SNP-only alignment)
   ├─ SNP-dists (Pairwise distance matrix)
   └─ IQ-TREE (ML phylogeny: GTR+G, 1000 ultrafast bootstrap)
      └─> 20,043 SNP positions used

8. Output & Results
   └─> Comprehensive phylogenetic trees, ANI matrices, AMR profiles
       sal(E) genetic context documentation
```

---

## Requirements

### Software
All tools are managed via conda/mamba environments:

```bash
mamba create -n download        -c bioconda ncbi-datasets-cli
mamba create -n checkm2         -c bioconda checkm2
mamba create -n quast           -c bioconda quast=5.2.0
mamba create -n bakta           -c bioconda bakta=1.12.0
mamba create -n fastani         -c bioconda fastani
mamba create -n panaroo         -c bioconda panaroo
mamba create -n gubbins         -c bioconda gubbins
mamba create -n clonalframe     -c bioconda clonalframeml
mamba create -n snpsites        -c bioconda snp-sites snp-dists
mamba create -n iqtree          -c bioconda iqtree
mamba create -n amrfinder       -c bioconda ncbi-amrfinderplus
mamba create -n card_new        -c bioconda rgi
mamba create -n abricate        -c bioconda abricate
mamba create -n virulencefinder -c bioconda virulencefinder
mamba create -n alignment       -c bioconda mafft iqtree
mamba create -n rotifer         -c bioconda rotifer
```

### Databases
| Database | Version | Path |
|---|---|---|
| Bakta full DB | v6.0 (2025-02-24) | `/databases/bakta_db/db` |
| CheckM2 | 2026-05-15.1 | `/databases/CheckM2_database/uniref100.KO.1.dmnd` |
| AMRFinderPlus | 2026-05-15.1 | auto-managed |
| CARD | 3.2.7 | `/databases/card/card.json` |
| VirulenceFinder DB | - | `/databases/virulencefinder_db` |
| Pfam | v35.0 | auto-managed / `/databases/pfam/Pfam-A.hmm` |

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

jq -r '.accession + "\t" + (.assemblyInfo.biosample.accession // "NO_BIOSAMPLE")' \
    assembly_data_report.jsonl > /tmp/assembly_biosample.txt

declare -A best_acc best_ver
while IFS=$'\t' read -r accession biosample; do
    # see scripts/dedup.sh for full logic
done < /tmp/assembly_biosample.txt
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
conda activate quast_5.2.0
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

# All vs reference genome (GCF_002442935.1 = NCBI reference strain)
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

# Input: Gubbins recombination-filtered SNP alignment (20,043 SNPs)
# Constant sites correction applied with -fconst
iqtree \
    -s gubbins_results/nepalensis.filtered_polymorphic_sites.fasta \
    -m GTR+G \
    -fconst 590159,246851,316673,503660 \
    -bb 1000 \
    -nt 20 \
    -pre iqtree_results/nepalensis

# Results:
#   Input:               Gubbins recombination-filtered SNPs
#   SNP positions:       20,043
#   Constant sites:      1,657,343 (corrected via -fconst)
#   Bootstrap:           1,000 UFBoot
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
# Identical genome pairs detected:
#   GCF_900458695.1 = GCF_014635045.1
#   GCF_051590305.1 = GCF_051590285.1
```

### Step 15 — Antimicrobial Resistance Analysis

#### AMRFinderPlus & CARD
```bash
conda activate amrfinder

mkdir -p amrfinder_results
for faa in bakta_results/*/*.faa; do
    accession=$(basename "$faa" .faa)
    [[ "$accession" == *.hypotheticals ]] && continue
    amrfinder \
        -p "$faa" \
        --organism Staphylococcus_aureus \
        --output amrfinder_results/${accession}.tsv \
        --threads 4
done

# Core resistome (100% genomes): sal(E)
# Accessory: fosB (64%), tet(K) (25%), mph(C) (17%), mecA (6%)
```

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
```

#### sal(E) Genetic Context Analysis

##### Extract sal(E) Flanking Sequences (ROTIFER)
```bash
conda activate rotifer

mkdir -p gene_extraction/sale_context

# Use ROTIFER to extract sal(E) and flanking sequences from all 36 genomes
# ROTIFER (Rapid Open-source Tools and Infrastructure for Data Exploration and Research)
# Reference: https://github.com/leepbioinfo/rotifer

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

# Translate flanking sequences (if needed) or use nucleotide hmmscan
# against Pfam database to identify conserved domains in vicinity of sal(E)

hmmscan \
    --domtblout gene_extraction/sale_context/pfam_annotations/flanking.domtbl \
    /databases/pfam/Pfam-A.hmm \
    gene_extraction/sale_context/*_flanking.faa

# Alternative: Extract ORFs from flanking regions, translate, and annotate
python3 extract_orfs_flanking.py \
    gene_extraction/sale_context/ \
    gene_extraction/sale_context/pfam_annotations/

# Results:
#   - Identification of genetic organization (operons, adjacent genes)
#   - Pfam domain annotations for neighboring coding sequences
#   - Assessment of sal(E) genomic context conservation
```

---

### Step 16 — Virulence Factor Detection

#### VirulenceFinder
```bash
conda activate virulencefinder

# Download and index database
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
python3 extract_sale.py
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
#   Conclusion: sal(E) is highly conserved, single ancestral acquisition
```

### Step 18 — Merge All Results (Python)
```python
# Scripts available:
# merge_fastani.py         -- adds FastANI vs reference to high-quality table
# merge_amr_tables_v2.py   -- adds AMRFinder + CARD results
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
scp -r zbecerra@davinci.icb.usp.br:/path/to/quast_results \
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
│   ├── nepalensis.filtered_polymorphic_sites.fasta  # input for SNP-sites + IQ-TREE
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
│   ├── nepalensis.treefile       # ML tree based on Gubbins SNPs
│   ├── nepalensis.contree        # bootstrap consensus tree
│   ├── nepalensis.iqtree         # full report
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
│   ├── sale_context/                        # ROTIFER extracted context
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

| Analysis | Result |
|---|---|
| Total genomes downloaded | 88 |
| Unique genomes after deduplication | 46 |
| High-quality genomes | 36 |
| Core genome size | 2,348 genes |
| Soft core genes | 40 genes |
| Shell genes | 485 genes |
| Cloud genes | 1,880 genes |
| Total pangenome | 4,753 genes |
| Pangenome type | Open |
| Gubbins SNPs (total) | 32,454 |
| SNPs used for phylogeny (IQ-TREE) | 20,043 |
| Recombination blocks (Gubbins) | 211 |
| SNPs from recombination | 1,908 (5.88%) |
| Gubbins r/m | 0.100 |
| ClonalFrameML R/theta | 3.97 x 10-5 |
| ClonalFrameML delta | 97.58 bp |
| ClonalFrameML nu | 2.09 x 10-4 |
| ClonalFrameML r/m | 8.09 x 10-7 |
| Core resistome | sal(E) |
| Accessory resistance genes | fosB, tet(K), mph(C), mecA |
| sal(E) conservation | 97.8% (531/543 sites) |
| sal(E) genetic context | Highly conserved across all genomes |
| Virulence genes detected | None |

---

## sal(E) Genetic Context Analysis

The sal(E) gene and its flanking sequences were extracted from all 36 genomes using
ROTIFER (Rapid Open-source Tools and Infrastructure for Data Exploration and Research;
https://github.com/leepbioinfo/rotifer). Approximately 2 kb regions flanking sal(E) on
both sides were retrieved to determine the genetic organization and evolutionary
conservation of the locus.

The extracted flanking sequences and their associated nucleotide sequences were annotated
using the Pfam database to identify conserved protein domains in adjacent coding regions.
This analysis revealed:

- **Genetic organization**: Consistent synteny across all S. nepalensis isolates
- **Associated genes**: Identification of neighboring genes and their functional domains
- **Evolutionary conservation**: Evidence of single ancestral acquisition with minimal
  subsequent modification
- **Functional context**: Coordination of sal(E) with neighboring metabolic or regulatory
  genes

---

## Recombination Parameters

```
r/m = R/theta x delta x nu
r/m = 3.97e-05 x 97.58 x 2.09e-04
r/m = 8.09e-07

Interpretation: S. nepalensis is a predominantly clonal species with
very low recombination rates, consistent with other CoNS species.
Short recombinant tracts (97.58 bp) and low divergence suggest
within-species recombination rather than external imports.
```

---

## Phylogenetic Analysis Notes

The maximum likelihood phylogeny was inferred using IQ-TREE v3.1.2 with the GTR+G
substitution model. The input alignment consisted of recombination-filtered SNPs
identified by Gubbins (20,043 SNPs from 36 genomes). Constant sites
(A=590,159; C=246,851; G=316,673; T=503,660) were included for branch length
correction using the -fconst option. Bootstrap support was assessed with
1,000 ultrafast bootstrap replicates.

---

## Virulence Factor Analysis

No virulence genes were detected in any of the 36 S. nepalensis genomes using:
- VirulenceFinder (databases: s.aureus_toxin, s.aureus_exoenzyme, s.aureus_hostimm;
  identity >= 90%, coverage >= 60%)
- Abricate with VFDB (identity >= 80%, coverage >= 60%)

This confirms the non-pathogenic nature of S. nepalensis and its distinction
from pathogenic staphylococci such as S. aureus.

---

## Software Versions

| Tool | Version | Reference |
|---|---|---|
| NCBI Datasets CLI | 18.29.1 | NCBI |
| CheckM2 | 1.1.0 | Chklovski et al. 2023 |
| QUAST | 5.2.0 | Gurevich et al. 2013 |
| Bakta | 1.12.0 | Schwengers et al. 2021 |
| FastANI | - | Jain et al. 2018 |
| Panaroo | - | Tonkin-Hill et al. 2020 |
| Gubbins | - | Croucher et al. 2015 |
| ClonalFrameML | 1.20 | Didelot & Wilson 2015 |
| SNP-sites | 2.5.1 | Page et al. 2016 |
| SNP-dists | 1.2.0 | Seemann 2018 |
| IQ-TREE | 3.1.2 | Minh et al. 2020 |
| AMRFinderPlus | 4.2.7 | Feldgarden et al. 2021 |
| RGI/CARD | - (DB: 3.2.7) | Alcock et al. 2023 |
| ROTIFER | - | https://github.com/leepbioinfo/rotifer |
| Pfam | v35.0 | Mistry et al. 2021 |
| VirulenceFinder | - | Joensen et al. 2014 |
| Abricate | - | Seemann 2020 |
| MAFFT | 7.525 | Katoh & Standley 2013 |

---

## References

- Alcock et al. (2023) CARD 2023. Nucleic Acids Research
- Chklovski et al. (2023) CheckM2. Nature Methods
- Croucher et al. (2015) Gubbins. Nucleic Acids Research
- Didelot & Wilson (2015) ClonalFrameML. PLoS Computational Biology
- Feldgarden et al. (2021) AMRFinderPlus. Scientific Reports
- Gurevich et al. (2013) QUAST. Bioinformatics
- Jain et al. (2018) FastANI. Nature Communications
- Joensen et al. (2014) VirulenceFinder. Journal of Clinical Microbiology
- Katoh & Standley (2013) MAFFT. Molecular Biology and Evolution
- Minh et al. (2020) IQ-TREE 2. Molecular Biology and Evolution
- Mistry et al. (2021) Pfam. Nucleic Acids Research
- Page et al. (2016) SNP-sites. Microbial Genomics
- Schwengers et al. (2021) Bakta. Microbial Genomics
- Tonkin-Hill et al. (2020) Panaroo. Genome Biology

---

## Citation

If you use this pipeline please cite:

> Becerra et al. (2026). Phylogeographic analysis of Staphylococcus nepalensis reveals global occurrence of antimicrobial-resistant lineages carrying the intrinsic sal(E) resistance gene.

---

## License
MIT License
