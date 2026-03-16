# Normalization Techniques in Omics Data Analysis
### BIOINF 504: Rigor & Reproducibility in Bioinformatics
##### Created with the assistance of claude.ai
---

> **Learning Objectives:** By the end of this module, students should be able to (1) explain why normalization is necessary in omics data, (2) select an appropriate normalization strategy for a given data type and experimental design, (3) implement normalization using standard R/Python packages, and (4) critically evaluate the pros and cons of each approach.

---

## Table of Contents

1. [Why Normalization Matters](#why-normalization-matters)
2. [Genomics: DNA-Seq & ChIP-Seq Normalization](#1-genomics-dna-seq--chip-seq-normalization)
3. [Transcriptomics: Bulk RNA-Seq Normalization](#2-transcriptomics-bulk-rna-seq-normalization)
4. [Transcriptomics: Single-Cell RNA-Seq Normalization](#3-transcriptomics-single-cell-rna-seq-scrna-seq-normalization)
5. [Proteomics Normalization](#4-proteomics-normalization)
6. [Metabolomics Normalization](#5-metabolomics-normalization)
7. [Choosing the Right Method: A Decision Guide](#choosing-the-right-method-a-decision-guide)
8. [Student Exercise: scRNA-seq Normalization in Practice](#student-exercise-scrna-seq-normalization-in-practice)

---

## Why Normalization Matters

Raw omics measurements are confounded by **technical variation** that has nothing to do with biology. Sources of technical noise include:

- **Sequencing depth** (library size): A sample sequenced to 50M reads will have ~2× higher raw counts than the same sample at 25M reads, even if the biology is identical.
- **RNA composition bias**: If a few highly expressed genes dominate one sample, they consume read budget and artificially deflate counts of other genes.
- **Batch effects**: Samples processed on different days, by different technicians, or on different instruments systematically differ.
- **Cell-to-cell technical variation** (scRNA-seq): Droplet capture efficiency, ambient RNA contamination, and doublets all add noise.
- **Sample loading variation** (proteomics/metabolomics): Differences in total protein or metabolite input per sample.

Normalization aims to remove these technical biases while **preserving true biological signal**. Choosing the wrong method — or skipping normalization entirely — is one of the most common sources of irreproducible results in omics analyses.

---

## 1. Genomics: DNA-Seq & ChIP-Seq Normalization

In genomics, normalization is most relevant for **copy number variation (CNV) analysis**, **ChIP-seq peak calling**, and **ATAC-seq accessibility comparisons** where read depth and mappability differ across samples or genomic regions.

---

### 1.1 RPKM and FPKM (Reads/Fragments Per Kilobase per Million mapped reads)

**RPKM** (Mortazavi et al., *Nature Methods*, 2008) was introduced for **single-end** RNA-seq, where each sequenced read corresponds to one fragment. **FPKM** is the paired-end equivalent: in paired-end sequencing, two reads are sequenced from the two ends of a single cDNA fragment. FPKM counts the number of *fragments* (insert molecules), not individual reads, to avoid counting the same molecule twice. In single-end experiments, RPKM and FPKM produce identical values. The formulas are otherwise identical.

**Formula (RPKM):**

```
RPKM = (Read counts × 10^6) / (Feature length in kb × Total mapped reads)

Equivalently: RPKM = (Read counts) / (Feature length in kb × Total mapped reads in millions)
```

**Step-by-step order of operations for RPKM:**

1. Divide total mapped reads by 1,000,000 → per-million scaling factor
2. Divide raw read counts by this scaling factor → Reads Per Million (RPM); corrects for sequencing depth
3. Divide RPM values by gene length in kilobases → RPKM; corrects for gene length

**Conceptual Example:**

Gene A is 2,000 bp long and has 1,000 reads in a sample with 10 million mapped reads.

```
RPKM = (1,000 × 10^6) / (2,000 bp × 10,000,000 reads)
     = 1,000 / (2.0 kb × 10 M)
     = 50
```

Gene B is 500 bp long and has 1,000 reads in the same sample.

```
RPKM = 1,000 / (0.5 kb × 10 M) = 200
```

Even though both genes have identical raw counts, Gene B appears more highly expressed *per unit length*, correcting for the fact that longer genes accumulate more reads by chance (longer transcripts have more positions for reads to land on).

**Why RPKM and FPKM are not reliably comparable across samples:**

The critical flaw of RPKM/FPKM is that the *sum* of all RPKM (or FPKM) values differs from sample to sample. This is a direct consequence of the order of operations: depth normalization is performed *before* length normalization. As a result, RPKM values are internally inconsistent as a measure of relative molar concentration: if gene A has RPKM = 50 in Sample 1 and RPKM = 50 in Sample 2, you cannot conclude that the same proportion of transcripts came from gene A in both samples — the denominators used to compute each sample's values are not equivalent (Wagner, Kin & Lynch, *Theory in Biosciences*, 2012).

**Pros:**
- Simple and relatively interpretable
- Corrects for both sequencing depth and feature length within a sample
- Suitable for comparing genes within the same sample

**Cons:**
- **Not reliably comparable across samples**: the sum of RPKM/FPKM values differs between samples, so a given RPKM value does not represent the same relative abundance in two different libraries (Wagner et al., 2012)
- Sensitive to composition bias from a small number of highly expressed genes
- RPKM/FPKM ≠ mRNA molar concentration; they are proportional to relative mass (nucleotides), not relative molarity

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `edgeR` | `rpkm(dge_object, gene.length = gene_lengths)` |
| R | `DESeq2` | Not built-in; typically computed manually after `counts()` |
| Python | `pydeseq2` | Manual: `counts / (lengths_kb * total_counts_M)` |

---

### 1.2 Coverage Normalization (Reads Per Million — RPM / CPM)

**Formula:**

```
CPM = (Raw counts × 10^6) / Total mapped reads in sample
```

Used for ChIP-seq, ATAC-seq, and CUT&RUN where **feature length differences are less relevant** (all peaks are roughly similar in size, or you are comparing signal at fixed genomic windows).

**Example:**

A ChIP-seq peak has 500 reads in a sample with 20M total reads.

```
CPM = (500 × 10^6) / 20,000,000 = 25 CPM
```

The same peak in a deeper sample (40M reads) with 1,000 reads:

```
CPM = (1,000 × 10^6) / 40,000,000 = 25 CPM ✓ (comparable)
```

**Pros:**
- Very simple; fast to compute
- Makes samples with different sequencing depth comparable
- Standard for ChIP-seq signal track generation (bigWig files)

**Cons:**
- Does not correct for feature/region length
- Susceptible to composition bias (if a few peaks dominate)
- Not appropriate for differential expression without further modeling

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `edgeR` | `cpm(dge_object, log = FALSE)` |
| R | `DESeq2` | Manual: `counts(dds) / colSums(counts(dds)) * 1e6` |
| Python | `deeptools` (CLI) | `bamCoverage --normalizeUsing CPM` |
| Python | `pybedtools` | Custom implementation over coverage arrays |

---

## 2. Transcriptomics: Bulk RNA-Seq Normalization

Bulk RNA-seq measures averaged gene expression across thousands to millions of cells. The main normalization goals are: (1) correct for sequencing depth, and (2) correct for composition bias so that differential expression analysis is valid.

---

### 2.1 TPM (Transcripts Per Million)

TPM was proposed by Li & Dewey (2011) and formalized by Wagner, Kin & Lynch (2012) as a correction to the cross-sample inconsistency of RPKM/FPKM. The fix is elegant: **reverse the order of operations** — normalize for gene length *first*, then normalize for sequencing depth *second*.

**Formula (two-step):**

```
Step 1: RPK  = Read counts / Feature length in kb        ← length normalization first
Step 2: TPM  = (RPK / Sum of all RPKs in sample) × 10^6  ← then depth normalization
```

**Why the order matters — the key difference from RPKM/FPKM:**

In RPKM, depth is corrected first, then length. This means the per-million scaling factor is computed from *raw* read counts, and the resulting RPKM values for all genes do not sum to a fixed constant — their sum varies from sample to sample.

In TPM, length is corrected first (producing RPK values). The *sum of all RPK values* becomes the per-million denominator. Because every gene's length-corrected count is divided by the same within-sample total, all TPM values in a given sample sum to exactly 10^6 — and this holds for every sample equally.

**This fixed sum is what makes TPM cross-sample comparable:** if Gene A has TPM = 50 in Sample 1 and TPM = 50 in Sample 2, you know that Gene A contributed the same proportion of the total length-corrected read pool in both samples. With RPKM, this interpretation is not valid because the denominators differ between samples.

**Example:**

| Gene | Counts | Length (kb) | RPK (= Counts / Length)  |
|------|--------|-------------|--------------------------|
| A    | 1,000  | 2.0         | 500                      |
| B    | 500    | 0.5         | 1,000                    |
| **Sum RPK** | | | **1,500**           |

```
TPM_A = (500 / 1,500) × 10^6 ≈ 333,333
TPM_B = (1,000 / 1,500) × 10^6 ≈ 666,667
Sum   = 1,000,000 ✓  (always exactly 10^6 within any sample)
```

For comparison, if these were computed as RPKM in a 10M-read library:

```
RPKM_A = 1,000 / (2.0 × 10) = 50
RPKM_B = 500  / (0.5 × 10)  = 100
Sum     = 150  (varies between samples — not a fixed constant)
```

> **Note on cross-sample comparability:** While TPM is better behaved than RPKM/FPKM, it is still a relative measure. Both RPKM/FPKM and TPM represent the abundance of a transcript *relative to all other sequenced transcripts in that sample*, and therefore depend on the composition of the RNA population. If two samples differ strongly in their transcriptome composition (e.g., different cell types, or different library preparation protocols), direct TPM comparisons can still be misleading (Zhao, Ye & Stanton, *RNA*, 2020). For rigorous differential expression analysis, count-based methods (DESeq2, edgeR) are preferred.

**Pros:**
- The sum of all TPM values is always 10^6 per sample — this fixed denominator makes within-metric proportions directly comparable across samples (unlike RPKM/FPKM)
- Corrects for both sequencing depth and gene length within a sample
- Preferred for reporting gene expression levels and for use with transcript quantification tools (Salmon, kallisto)

**Cons:**
- Still a relative, compositional measure — affected by the RNA composition of the sample
- Not suitable as direct input to DESeq2 or edgeR (these require raw integer counts to model dispersion correctly)
- Transcript-level quantification requires careful length calculation; use `tximport` to correctly import TPM from Salmon/kallisto into Bioconductor workflows

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `edgeR` | Manual from `rpkm()` output; normalize RPK by column sums × 1e6 |
| R | `tximport` | `tximport(files, type="salmon", countsFromAbundance="scaledTPM")` |
| Python | `pyroe` / manual | Load directly from Salmon/kallisto output, or compute manually |

---

### 2.2 TMM (Trimmed Mean of M-values) — edgeR

TMM is a **between-sample normalization** method that assumes the majority of genes are NOT differentially expressed. It estimates a scaling factor for each sample relative to a reference sample by trimming the most extreme log-fold changes (M-values) and extreme average expression values (A-values), then computing a **weighted** trimmed mean of the remaining M-values (weights are inverse-variance, computed via the delta method for log-binomial data). By default, `edgeR`'s `calcNormFactors` trims the top and bottom **30%** of M-values and the top and bottom **5%** of A-values (Robinson & Oshlack, *Genome Biology*, 2010; confirmed in Maza et al., *Frontiers in Genetics*, 2016).

**Conceptual Example:**

Imagine two RNA-seq samples where, in Sample 2, a few highly induced genes (e.g., viral transcripts) dominate the library and consume a disproportionate fraction of reads. Without correction, all other genes appear "down-regulated" in Sample 2 — but this is a composition artifact. TMM corrects for this:

1. A reference sample is chosen (the library whose upper-quartile count-per-million is closest to the mean upper-quartile across all libraries).
2. For each gene, compute **M-values** (log₂ fold-change between the test and reference) and **A-values** (mean log₂ expression across both).
3. Trim the top and bottom 30% of M-values and top/bottom 5% of A-values to remove extreme genes.
4. Compute the inverse-variance-weighted mean of the remaining M-values. This is the log₂ normalization factor for that sample.
5. The resulting **normalization factors** scale the effective library sizes, not the counts directly.

**Pros:**
- Robust to outlier genes (trimming removes compositional extremes)
- Default and well-validated method in `edgeR`
- Works well when <50% of genes are differentially expressed

**Cons:**
- Assumes majority of genes are not DE — violated in very distinct cell types or extreme perturbations
- Produces scaling factors, not transformed values; downstream tools must use these factors
- Less intuitive than TPM

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `edgeR` | `calcNormFactors(dge_object, method = "TMM")` |
| R | `edgeR` | Access factors via `dge$samples$norm.factors` |
| Python | `pydeseq2` | Not directly implemented; use R via `rpy2` |

---

### 2.3 DESeq2 Median-of-Ratios Normalization

DESeq2 computes a **geometric mean** for each gene across all samples (serving as a reference), then divides each sample's gene counts by this reference to get per-gene ratios. The **median** of all these ratios becomes the size factor for that sample.

**Step-by-step for 2 genes, 2 samples:**

```
           Gene1  Gene2
Sample A:    100    200
Sample B:    200    400

Geometric mean: Gene1 = sqrt(100*200) ≈ 141.4
                Gene2 = sqrt(200*400) ≈ 282.8

Ratios:
  Sample A: Gene1 = 100/141.4 = 0.707,  Gene2 = 200/282.8 = 0.707 → median = 0.707
  Sample B: Gene1 = 200/141.4 = 1.414,  Gene2 = 400/282.8 = 1.414 → median = 1.414

Normalized counts = Raw counts / size factor
  Sample A (÷ 0.707): Gene1 = 141, Gene2 = 283
  Sample B (÷ 1.414): Gene1 = 141, Gene2 = 283 ✓
```

**Pros:**
- Highly robust to extreme outlier genes (median is used, not mean)
- Handles RNA composition bias effectively
- Integrated into the DESeq2 negative binomial model
- Well-validated in hundreds of publications

**Cons:**
- Requires count data (integers); cannot be applied to pre-normalized values
- Genes with zero counts in any sample are excluded from size factor estimation (problematic for sparse data like scRNA-seq)
- Not designed for cross-experiment comparisons

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `DESeq2` | `estimateSizeFactors(dds)` then `counts(dds, normalized = TRUE)` |
| Python | `pydeseq2` | `dds.fit_size_factors()` then access via `dds.layers["normed_counts"]` or `dds.obsm["size_factors"]` |

---

### 2.4 Quantile Normalization

Forces all samples to have the same distribution by replacing values with the average value of the same rank across samples. Originally developed for microarrays.

**Simple Example (3 genes, 2 samples):**

```
Raw:      S1=[10, 30, 20]   S2=[50, 20, 10]
Sorted:   S1=[10, 20, 30]   S2=[10, 20, 50]
Means:    rank1=(10+10)/2=10, rank2=(20+20)/2=20, rank3=(30+50)/2=40
Assign means back by original rank:
  S1: 10→10, 30→40, 20→20  → [10, 40, 20]
  S2: 50→40, 20→20, 10→10  → [40, 20, 10]
```

**Pros:**
- Equalizes distributions across all samples
- Commonly used for microarray and methylation array data
- Implemented in many popular packages

**Cons:**
- **Aggressive**: assumes all samples should have the same distribution — biologically unrealistic for RNA-seq
- Can introduce artificial correlations between samples
- Generally **not recommended for RNA-seq** (distorts true abundance differences)
- Inappropriate when true global expression differences exist

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `preprocessCore` (Bioconductor) | `normalize.quantiles(matrix)` |
| R | `limma` | `normalizeQuantiles(matrix)` |
| Python | `sklearn.preprocessing` | Manual rank-based implementation |
| Python | `qnorm` | `qnorm.quantile_normalize(df)` |

---

## 3. Transcriptomics: Single-Cell RNA-Seq (scRNA-seq) Normalization

scRNA-seq presents **unique normalization challenges** absent in bulk RNA-seq:

- **Extreme sparsity**: Most gene-cell entries are zero (dropout events)
- **High technical variability**: Capture efficiency varies across cells
- **Cell-type heterogeneity**: Cells differ in total RNA content biologically
- **No "housekeeping" genes**: Expression is highly variable even for stable genes

---

### 3.1 Library Size (Total Count) Normalization — "LogNormalize"

The simplest and most widely used scRNA-seq normalization. Each cell's counts are divided by its total count (library size), multiplied by a scale factor (typically 10,000 = "counts per 10K", or CP10K), then log-transformed.

**Formula:**

```
Normalized_ij = log( (count_ij / total_counts_j) × scale_factor + 1 )
```

**Example:**

Cell A has total counts = 5,000; Gene X count = 50
Cell B has total counts = 15,000; Gene X count = 150

```
Cell A norm = log((50/5000) × 10000 + 1) = log(100 + 1) ≈ 4.62
Cell B norm = log((150/15000) × 10000 + 1) = log(100 + 1) ≈ 4.62 ✓
```

Both cells show equal expression of Gene X after normalization — correct, since 50/5000 = 150/15000 = 1% of transcriptome.

**Pros:**
- Extremely simple and fast
- Default in Seurat; well-understood behavior
- log1p transformation compresses the dynamic range and partially stabilizes variance, making data more amenable to linear methods like PCA

**Cons:**
- Assumes all cells have the same total RNA content — biologically questionable (large cells, proliferating cells genuinely have more RNA)
- Sensitive to "doublets" (two cells captured together with inflated library sizes)
- Does not fully correct for zero-inflation or technical dropout

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `Seurat` | `NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)` |
| Python | `scanpy` | `sc.pp.normalize_total(adata, target_sum=1e4)` then `sc.pp.log1p(adata)` |

---

### 3.2 Scran Pooling-Based Normalization

Developed specifically for scRNA-seq, `scran` addresses the problem that **per-cell size factor estimation fails** for sparse data (many zeros make the DESeq2 median-of-ratios approach unreliable at the single-cell level). 

The key idea: **pool cells into groups, normalize the pool (reducing sparsity), then deconvolve** pool-based size factors back to individual cells using a system of linear equations.

**Conceptual Steps (Lun, Bach & Marioni, *Genome Biology*, 2016):**

1. Pre-cluster cells roughly (e.g., with `quickCluster`) so that pools are drawn from cells with broadly similar expression profiles, reducing the number of DE genes within any pool.
2. Within each cluster, cells are **sorted by increasing library size**. A **sliding window** of varying sizes (by default, pool sizes ranging from ~21 to ~101 cells at intervals of 5) is applied to this ordering.
3. At each window position, expression counts are **summed** across the pooled cells, reducing sparsity and enabling reliable median-of-ratios normalization on the pool.
4. A pool-based size factor is estimated as the **median ratio** between the pooled counts and a reference pseudo-cell (the average of all cells). Because the pool's size factor equals the sum of the constituent cells' size factors, each pool defines one linear equation in the unknown cell-level size factors.
5. Pooling across many window positions yields an **over-determined linear system**, solved by least-squares to recover per-cell size factors. Size factors computed within each cluster are then rescaled to be comparable across clusters.

**Pros:**
- More accurate size factors than simple library size for heterogeneous populations
- Handles zeros better than DESeq2 median-of-ratios applied naively to single cells
- Accounts for cluster-specific expression profiles
- Well-validated and widely recommended for UMI-based scRNA-seq datasets (e.g., in the Bioconductor OSCA workflow)

**Cons:**
- More complex and computationally intensive than `LogNormalize`
- Requires pre-clustering (introduces a circularity: normalization depends on clustering, clustering depends on normalization)
- Assumes size factors are cell-specific (not gene-specific)
- Less intuitive to explain and debug

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `scran` (Bioconductor) | `computeSumFactors(sce, clusters = pre.clusters)` |
| R | `scran` | Pre-clustering: `quickCluster(sce)` |
| Python | `scanpy` + `rpy2` or `anndata2ri` | Call R `scran` via `rpy2`/`anndata2ri` bridge; native Python port `scranPY` also exists |

---

### 3.3 SCTransform (Regularized Negative Binomial Regression)

SCTransform (Hafemeister & Satija, 2019) takes a fundamentally different approach: rather than dividing by library size, it **regresses out technical variation** using a regularized negative binomial regression for each gene.

**Model:** For each gene `g` in cell `j`:

```
counts_gj ~ NB(mu_gj, theta_g)
log(mu_gj) = beta0_g + beta1_g × log(library_size_j) + [other covariates]
```

The **Pearson residuals** from this model are used as normalized values — they represent how much each gene deviates from what would be expected given the cell's library size.

**Regularization:** Because fitting separate NB models for ~20,000 genes is noisy, SCTransform pools information across genes with similar expression levels (kernel regression over log-mean expression) to get stable parameter estimates.

**SCTransform v2** (2022) additionally clips residuals and uses 5,000 variable genes selected during the GLM fit itself.

**Pros:**
- Does not assume equal total RNA per cell
- Naturally handles the mean-variance relationship in scRNA-seq
- Can regress out known confounders (cell cycle, mitochondrial %) simultaneously
- Often produces cleaner UMAP/clustering results

**Cons:**
- Residuals are not counts — cannot be used with count-based DE methods (use corrected counts for DE)
- Computationally slower than LogNormalize (though parallelizable)
- Pearson residuals are not easily interpretable as "expression levels"
- Can overcorrect in datasets with very strong biological signals

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `Seurat` (v4+) | `SCTransform(seurat_obj, vst.flavor = "v2", vars.to.regress = "percent.mt")` |
| Python | `scanpy` | Not natively implemented; use `rpy2`/`anndata2ri` to call R's `Seurat::SCTransform`, or use `scanpy`'s `sc.experimental.pp.highly_variable_genes` with `flavor="seurat_v3"` as an approximation |

---

### 3.4 Analytic Pearson Residuals (glmGamPoi)

A fast analytical approximation of SCTransform residuals, avoiding the computationally expensive iterative fitting. Uses a generalized linear model with gamma-Poisson (NB) distribution.

**Pros:**
- Much faster than SCTransform (10–50×)
- Good approximation of SCTransform residuals
- Integrated into current Seurat v5 default workflow

**Cons:**
- Same conceptual limitations as SCTransform regarding residual interpretation

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `glmGamPoi` | `glm_gp(counts_matrix, design = ~log_umi)` |
| R | `Seurat` | `SCTransform(obj, vst.flavor="v2")` (uses `glmGamPoi` internally) |

---

## 4. Proteomics Normalization

Proteomics data (from mass spectrometry) has different properties than RNA-seq: intensities are continuous (not counts), many missing values exist (proteins not detected), and dynamic range is enormous (~4–6 orders of magnitude).

---

### 4.1 Total Ion Current (TIC) / Sum Normalization

Divide each sample's protein intensities by the total sum of all intensities in that sample (analogous to CPM in RNA-seq). Often applied to label-free quantification (LFQ) data.

**Example:**

```
Sample 1 intensities: [100, 200, 300, 400] → Sum = 1,000
Sample 2 intensities: [200, 400, 600, 800] → Sum = 2,000

Normalized (× 10^6 / sum):
Sample 1: [100, 200, 300, 400]
Sample 2: [100, 200, 300, 400] ✓ (identical biological composition)
```

**Pros:**
- Conceptually simple; widely applicable
- Standard for label-free proteomics
- Implemented in most proteomics software platforms

**Cons:**
- Highly sensitive to high-abundance "contaminant" proteins (keratins, albumin)
- Assumes total protein is constant across samples — inappropriate if samples have genuine differences in total protein
- Does not handle missing values

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `DEP` (Bioconductor) | `normalize_vsn(se_object)` for VSN; manual sum normalization via `sweep(assay(se), 2, colSums(assay(se)), "/")` |
| R | `MSnbase` | `normalise(mse, method = "sum")` |
| Python | `pandas` / manual | `df.div(df.sum(axis=0), axis=1)` (column-wise sum normalization) |

---

### 4.2 Median Centering / Median of Ratios (LFQ)

Shift each sample's log-intensities so that the **median log-ratio** to a reference (or grand median) is zero. Analogous to DESeq2 size factors.

**Example:**

```
Log2 intensities (after missing value imputation):
Protein   S1    S2    S3
A        10.0  11.5  10.8
B         8.0   9.5   8.8
C        12.0  13.5  12.8

Per-sample median: 10.0  11.5  10.8
Grand median: 10.77

Shift: S1 += +0.77, S2 += -0.73, S3 += -0.03
```

**Pros:**
- Robust to outlier proteins
- Simple and transparent
- MaxQuant's LFQ algorithm uses a variant of this approach

**Cons:**
- Requires imputation of missing values first (which introduces its own biases)
- Assumes majority of proteins are not differentially abundant (same caveat as TMM)

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `MSnbase` | `normalise(mse, method = "center.median")` |
| R | `DEP` (Bioconductor) | `normalize_vsn(se)` implements VSN (see §4.3); for median centering use `normalise(mse, method="center.median")` via `MSnbase` |
| Python | `pandas` / manual | `df.subtract(df.median(axis=0), axis=1)` on log-scale |

---

### 4.3 VSN (Variance Stabilizing Normalization)

VSN (Huber et al., *Bioinformatics*, 2002) models the relationship between mean intensity and variance using an **arsinh (arcsinh) transformation** applied after an affine transformation of the raw intensities. The parameters are estimated from the data by maximum likelihood to make variance approximately constant across the dynamic range (variance stabilization).

**Transformation:**

```
h(x) = arcsinh(a + b × x)
```

where `a` (offset) and `b` (scale) are fit per sample by maximum likelihood to minimize the dependence of variance on mean intensity. The transformation simultaneously normalizes across samples and stabilizes variance — unlike log transformation, which only stabilizes variance for high-intensity features but not for low-intensity ones near the noise floor.

**Pros:**
- Simultaneously normalizes and stabilizes variance
- Works on raw intensities without log-transformation
- Handles low-intensity regions better than log-based methods

**Cons:**
- Less interpretable transformation
- Computationally requires iterative fitting
- Assumes a two-component variance model

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `vsn` (Bioconductor) | `justvsn(intensity_matrix)` |
| R | `DEP` | `normalize_vsn(se_object)` |

---

## 5. Metabolomics Normalization

Metabolomics data (from LC-MS or NMR) faces challenges including signal drift over long acquisition runs, ionization suppression, and extreme dynamic range. Multiple normalization layers are typically applied sequentially.

---

### 5.1 Internal Standard (IS) Normalization

**Spike in** one or more known stable-isotope-labeled compounds (e.g., deuterium-labeled metabolites) at a fixed concentration before sample preparation. Divide all metabolite intensities by the internal standard's measured intensity.

**Example:**

```
Deuterium-labeled glucose (d7-glucose) spiked at 100 ng/mL:

          Metabolite X    d7-glucose
Sample 1:    5,000          2,000
Sample 2:    8,000          4,000    ← higher instrument sensitivity this run

IS-normalized:
Sample 1: 5000 / 2000 = 2.5
Sample 2: 8000 / 4000 = 2.0  ← corrected for ionization variation
```

**Pros:**
- Directly corrects for extraction efficiency and ionization variation
- Most accurate method when appropriate IS compounds are available
- Handles run-to-run instrument drift

**Cons:**
- Expensive: requires synthesized isotope-labeled standards
- One IS may not represent all metabolites across different chemical classes
- Does not correct for sample-to-sample biological variation in loading

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `MetaboAnalystR` | `Normalization(mSet, "CompNorm", ...)` (compound normalization) |
| R | `xcms` | `PeakGroupsParam()` then manual IS division |
| Python | `pymzml` | Manual: `df.div(df['IS_column'], axis=0)` |

---

### 5.2 Total Signal / Probabilistic Quotient Normalization (PQN)

**Total Signal:** Divide each sample by the sum of all metabolite intensities. Analogous to TIC for proteomics.

**PQN** (Dieterle et al., *Analytical Chemistry*, 2006) was developed for NMR metabolomics but is applicable to LC-MS data. It is based on the observation that most metabolites are not changed between groups, and that the most likely explanation for differences between a sample and a reference is a single dilution-like factor.

**Steps:**
1. Compute a **reference spectrum**: the median intensity of each metabolite across all samples (after an initial median-fold scaling to put all samples on a comparable scale).
2. For each sample, compute the ratio of every metabolite's intensity to the corresponding reference value.
3. The **normalization factor** for that sample is the **median** of all these per-metabolite ratios — this is the "most probable quotient."
4. Divide each sample's intensities by its normalization factor.

This approach is more robust than simple total-signal normalization because the median ratio is insensitive to metabolites that are genuinely changed between conditions.

**Pros:**
- PQN is more robust than total-signal normalization
- Does not require internal standards
- Recommended for NMR metabolomics

**Cons:**
- Assumes majority of metabolites are not changed
- Sensitive to high-abundance metabolites

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `MetaboAnalystR` | `Normalization(mSet, rowNorm="ProbNorm", transNorm="NULL", scaleNorm="NULL", ref=NULL)` |
| R | `pmp` (Bioconductor) | `pqn_normalisation(data_matrix)` |
| Python | `pymetabo` / manual | Custom PQN implementation |

---

### 5.3 Quantile Normalization (Metabolomics)

Same approach as described for microarrays — forces all samples to the same distribution. More acceptable in metabolomics than RNA-seq when biological variation is expected to be small relative to technical variation.

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `preprocessCore` | `normalize.quantiles(matrix)` |
| Python | `qnorm` | `qnorm.quantile_normalize(df)` |

---

### 5.4 QC-Based Signal Drift Correction (LOESS / QCRSC)

In long metabolomics acquisition runs (dozens to hundreds of samples), **systematic signal drift** occurs as instrument sensitivity changes over time due to detector fouling, column degradation, or source contamination. To monitor and correct this, **pooled QC (PQC) samples** — typically a pooled mixture of all study samples — are injected at regular intervals throughout the run (e.g., every 5–10 study samples).

The correction strategy: fit a **LOESS (locally estimated scatterplot smoothing)** regression curve through the PQC signal intensities over injection order for each metabolite, then divide each study sample's intensity by the expected drift value at its injection position (interpolated from the fitted curve).

**Example:** If PQC injections show a steady 20% decline in a metabolite's signal from injection 1 to injection 100, LOESS fits this trend. A study sample injected at position 60 is corrected by dividing by the fitted LOESS value at position 60, removing the drift artifact.

**Pros:**
- Directly corrects for instrument drift — essential for large multi-day studies
- Non-parametric: no assumption about the functional form of drift
- PQC-based approach is widely adopted as best practice (e.g., MetaboLights, EMBL-EBI guidelines)

**Cons:**
- Requires regularly spaced PQC injections throughout the run (adds ~10–15% extra instrument time)
- LOESS smoothing parameter (span/`spar`) must be tuned; too small → overfitting noise; too large → undercorrection
- Cannot correct for between-sample biological loading differences (use PQN or IS normalization for that)
- Requires a minimum number of PQC replicates per metabolite for reliable fitting (typically ≥ 6–8 QC injections recommended)

**Packages & Functions:**

| Language | Package | Function |
|----------|---------|----------|
| R | `pmp` (Bioconductor) | `QCRSC(df, order, batch, classes, spar, minQC)` — QC-based LOESS signal correction |
| R | `limma` | `normalizeCyclicLoess(log_matrix)` — cyclic LOESS for between-sample normalization (not drift correction) |
| Python | `statsmodels` + manual | `statsmodels.nonparametric.smoothers_lowess.lowess(qc_signal, injection_order)` then interpolate to study samples |

---

## Choosing the Right Method: A Decision Guide

```
What is your data type?
│
├─ Bulk RNA-seq
│   ├─ Comparing expression within same experiment? → DESeq2 median-of-ratios or TMM
│   ├─ Reporting/visualization? → TPM
│   └─ Microarray / log-scale data? → Quantile or LOESS
│
├─ scRNA-seq
│   ├─ Quick exploratory analysis? → LogNormalize (Seurat default)
│   ├─ Accurate size factors for heterogeneous populations? → scran
│   ├─ Regress out confounders / best UMAP? → SCTransform v2
│   └─ Large dataset, need speed? → Analytic Pearson residuals (glmGamPoi)
│
├─ Proteomics (LC-MS)
│   ├─ LFQ data? → Median centering or VSN
│   ├─ TMT/iTRAQ? → Reporter-ion based (mean centering)
│   └─ Absolute quantification? → IS normalization
│
└─ Metabolomics
    ├─ Targeted panel, have IS? → IS normalization
    ├─ Untargeted, large cohort? → PQN + LOESS drift correction
    └─ NMR data? → PQN or total area normalization
```

---

## Student Exercise: scRNA-seq Normalization in Practice

### Overview

In this exercise, you will normalize a **simulated single-cell RNA-seq dataset** using three different methods — `LogNormalize`, `scran`, and `SCTransform` — and compare their effects on downstream analysis (highly variable gene selection, PCA, and UMAP). You will observe how normalization choice affects biological interpretation.

The simulated dataset contains **500 cells** from **3 biologically distinct cell types**, with realistic technical variation (variable library sizes, dropout, and a "problematic" cell type with genuinely higher total RNA content).

---

### Part A: The Simulated Dataset

The code below generates the dataset. **Do not modify this section** — it represents your "raw data."

```r
# ============================================================
# STUDENT EXERCISE: scRNA-seq Normalization Comparison
# BIOINF 504: Rigor & Reproducibility in Bioinformatics
# ============================================================

# ----- Install required packages (run once) -----
# install.packages("BiocManager")
# BiocManager::install(c("scran", "scater", "SingleCellExperiment"))
# install.packages(c("Seurat", "ggplot2", "patchwork", "dplyr", "tidyr"))

library(Seurat)
library(scran)
library(scater)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(42)  # For reproducibility — always set a seed!

# ============================================================
# PART A: SIMULATE A REALISTIC scRNA-seq DATASET
# ============================================================
# We simulate 500 cells from 3 cell types:
#   - Type 1 (n=200): "Neurons"     — moderate total RNA (~3,000 UMIs)
#   - Type 2 (n=200): "Astrocytes"  — moderate total RNA (~3,500 UMIs)
#   - Type 3 (n=100): "Microglia"   — HIGH total RNA (~8,000 UMIs)
#
# The high total RNA in microglia is BIOLOGICALLY REAL (they are phagocytic
# cells with high transcriptional activity), not a technical artifact.
# This creates an interesting challenge for library-size normalization.
#
# We simulate 500 genes: 50 are "marker genes" for each type,
# 350 are "housekeeping" genes expressed at similar levels across types.

n_cells_per_type <- c(200, 200, 100)  # Neurons, Astrocytes, Microglia
n_genes          <- 500               # Total genes
n_cells          <- sum(n_cells_per_type)

# ---- Simulate base expression profiles per gene per cell type ----
# Each cell type has a mean expression vector for all 500 genes
simulate_type_means <- function(n_genes, base_mean, marker_genes, marker_fold) {
  # base_mean: average expression for housekeeping genes (counts scale)
  means <- rnbinom(n_genes, mu = base_mean, size = 5)      # housekeeping genes
  means[marker_genes] <- means[marker_genes] * marker_fold # up-regulate markers
  means <- pmax(means, 0)  # no negative values
  return(means)
}

# Define marker gene indices for each cell type
neuron_markers    <- 1:50        # genes 1-50
astrocyte_markers <- 51:100      # genes 51-100
microglia_markers <- 101:150     # genes 101-150

# Simulate mean expression per gene for each cell type
neuron_means    <- simulate_type_means(n_genes, base_mean = 3,  neuron_markers,    fold = 8)
astrocyte_means <- simulate_type_means(n_genes, base_mean = 3,  astrocyte_markers, fold = 8)
microglia_means <- simulate_type_means(n_genes, base_mean = 7,  microglia_markers, fold = 8)
#                                                          ^                               
# NOTE: microglia base_mean = 7 vs 3 for others — this simulates
# genuinely higher total RNA content (bigger, more active cells)

# ---- Sample cells from each type's expression distribution ----
# Each cell's counts drawn from negative binomial (overdispersed Poisson)
# to mimic real scRNA-seq count distributions

sample_cells <- function(n_cells, gene_means, size_factor_range) {
  # size_factor_range: min/max of cell-specific capture efficiency
  # (technical variation in library size within a cell type)
  cell_sf <- runif(n_cells,
                   min = size_factor_range[1],
                   max = size_factor_range[2])  # per-cell technical size factors
  
  counts <- sapply(1:n_cells, function(i) {
    # Scale means by cell-specific technical factor, then draw from NB
    cell_means <- gene_means * cell_sf[i]
    rnbinom(length(gene_means), mu = cell_means, size = 2)
  })
  return(counts)  # rows = genes, cols = cells
}

# Sample counts for each cell type
neuron_counts    <- sample_cells(200, neuron_means,    size_factor_range = c(0.5, 1.5))
astrocyte_counts <- sample_cells(200, astrocyte_means, size_factor_range = c(0.6, 1.4))
microglia_counts <- sample_cells(100, microglia_means, size_factor_range = c(0.7, 1.3))

# ---- Assemble full count matrix ----
counts_matrix <- cbind(neuron_counts, astrocyte_counts, microglia_counts)

# Name rows and columns
rownames(counts_matrix) <- paste0("Gene", 1:n_genes)
colnames(counts_matrix) <- paste0("Cell", 1:n_cells)

# True cell type labels (in a real experiment, these are UNKNOWN)
true_labels <- c(rep("Neuron", 200),
                 rep("Astrocyte", 200),
                 rep("Microglia", 100))

# ---- Introduce dropout (technical zeros) ----
# In scRNA-seq, low-expressed transcripts are often not captured ("dropout")
# We zero out counts with probability inversely proportional to expression
dropout_prob <- function(count) {
  # Higher counts = lower probability of dropout
  # Uses a logistic model: P(dropout) = 1 / (1 + exp(count - 2))
  1 / (1 + exp(count - 2))
}

dropout_mask <- matrix(
  rbinom(n = prod(dim(counts_matrix)),
         size = 1,
         prob = 1 - dropout_prob(counts_matrix)),  # 1 = keep, 0 = dropout
  nrow = nrow(counts_matrix)
)
counts_matrix_dropout <- counts_matrix * dropout_mask

cat("Dataset summary:\n")
cat("  Cells:", ncol(counts_matrix_dropout), "\n")
cat("  Genes:", nrow(counts_matrix_dropout), "\n")
cat("  Sparsity (% zeros):",
    round(100 * mean(counts_matrix_dropout == 0), 1), "%\n")
cat("  Median library sizes by true type:\n")
cat("    Neuron:",    median(colSums(counts_matrix_dropout[, true_labels == "Neuron"])), "\n")
cat("    Astrocyte:", median(colSums(counts_matrix_dropout[, true_labels == "Astrocyte"])), "\n")
cat("    Microglia:", median(colSums(counts_matrix_dropout[, true_labels == "Microglia"])), "\n")
```

---

### Part B: Apply Three Normalization Methods

```r
# ============================================================
# PART B: APPLY THREE NORMALIZATION METHODS
# ============================================================

# ---- B1: Create Seurat object from raw counts ----
seurat_obj <- CreateSeuratObject(
  counts = counts_matrix_dropout,
  project = "NormalizationExercise",
  min.cells = 3,    # keep genes detected in ≥ 3 cells
  min.features = 50 # keep cells with ≥ 50 detected genes
)

# Add true labels as metadata (in real analysis, you wouldn't have these yet!)
seurat_obj$true_type <- true_labels[match(colnames(seurat_obj),
                                          paste0("Cell", 1:n_cells))]

cat("After QC filtering:", ncol(seurat_obj), "cells,",
    nrow(seurat_obj), "genes retained\n")

# ---- B2: Visualize library size distribution BEFORE normalization ----
# This is critical QC — always examine your data before normalizing!
lib_sizes <- colSums(GetAssayData(seurat_obj, layer = "counts"))

p_lib <- data.frame(
  LibSize  = lib_sizes,
  CellType = seurat_obj$true_type
) %>%
  ggplot(aes(x = CellType, y = LibSize, fill = CellType)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("Neuron" = "#4C72B0",
                                "Astrocyte" = "#DD8452",
                                "Microglia" = "#55A868")) +
  labs(title = "Library Sizes Before Normalization",
       subtitle = "Note the higher counts in Microglia — is this biology or artifact?",
       x = "Cell Type", y = "Total UMI Counts") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

print(p_lib)

# ============================================================
# METHOD 1: LogNormalize (Seurat default)
# Divide each cell by its total counts, multiply by 10,000, log1p-transform
# This ASSUMES all cells have the same total RNA content
# ============================================================

seurat_log <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000  # "counts per 10K" (CP10K) then log1p
)

# Select highly variable genes (HVGs) — genes that vary more than expected
# by chance; these will drive dimensionality reduction
seurat_log <- FindVariableFeatures(seurat_log,
                                    selection.method = "vst",
                                    nfeatures = 100)

# Scale data: zero-center and unit-variance each gene (for PCA)
seurat_log <- ScaleData(seurat_log)

# Run PCA
seurat_log <- RunPCA(seurat_log, npcs = 20, verbose = FALSE)

# Run UMAP for 2D visualization
seurat_log <- RunUMAP(seurat_log, dims = 1:15, verbose = FALSE)

# ============================================================
# METHOD 2: scran pooling-based normalization
# More sophisticated: pools cells to handle zeros, then deconvolves
# to per-cell size factors
# ============================================================

# scran uses SingleCellExperiment (SCE) objects, not Seurat
sce <- SingleCellExperiment(
  assays = list(counts = GetAssayData(seurat_obj, layer = "counts"))
)

# Step 1: Quick pre-clustering (needed by scran to handle cell-type differences)
# This is the "chicken-and-egg" problem: we cluster before normalizing!
pre_clusters <- quickCluster(sce, min.size = 50)
cat("scran pre-clusters identified:", nlevels(pre_clusters), "\n")

# Step 2: Estimate per-cell size factors using pooling-and-deconvolution
sce <- computeSumFactors(sce,
                          clusters = pre_clusters,
                          min.mean = 0.1)  # only use genes with mean ≥ 0.1

# Inspect the estimated size factors
cat("scran size factor summary:\n")
print(summary(sizeFactors(sce)))

# Step 3: Apply size factors and log-normalize
sce <- logNormCounts(sce)  # uses sizeFactors stored in the SCE object

# Transfer scran-normalized data back into a Seurat object for consistent downstream analysis
seurat_scran <- seurat_obj  # copy metadata, structure
# Replace normalized data with scran's log-normalized values
seurat_scran[["RNA"]]@layers[["data"]] <- logcounts(sce)[
  rownames(seurat_scran), colnames(seurat_scran)]

# Proceed with HVG selection, scaling, PCA, UMAP
seurat_scran <- FindVariableFeatures(seurat_scran,
                                      selection.method = "vst",
                                      nfeatures = 100)
seurat_scran <- ScaleData(seurat_scran)
seurat_scran <- RunPCA(seurat_scran, npcs = 20, verbose = FALSE)
seurat_scran <- RunUMAP(seurat_scran, dims = 1:15, verbose = FALSE)

# ============================================================
# METHOD 3: SCTransform v2
# Regresses out library size effect using regularized NB regression
# Outputs Pearson residuals as normalized values
# NOTE: SCTransform replaces NormalizeData + ScaleData + FindVariableFeatures
# ============================================================

seurat_sct <- SCTransform(
  seurat_obj,
  vst.flavor = "v2",         # use the improved 2022 version
  variable.features.n = 100, # select 100 HVGs (matching other methods)
  verbose = FALSE
)

# SCTransform automatically sets the default assay to "SCT"
# The Pearson residuals are in the "scale.data" slot
seurat_sct <- RunPCA(seurat_sct, npcs = 20, verbose = FALSE)
seurat_sct <- RunUMAP(seurat_sct, dims = 1:15, verbose = FALSE)
```

---

### Part C: Compare Results

```r
# ============================================================
# PART C: COMPARE NORMALIZATION METHODS
# ============================================================

# ---- C1: Compare size factors / effective normalization ----
# For LogNormalize, the "effective size factor" is just library size / 10000
log_sf    <- lib_sizes[colnames(seurat_obj)] / 10000
scran_sf  <- sizeFactors(sce)[colnames(seurat_obj)]

# Compare size factors by cell type
sf_df <- data.frame(
  CellID      = colnames(seurat_obj),
  CellType    = seurat_obj$true_type,
  LogNorm_SF  = log_sf,
  Scran_SF    = scran_sf
)

# Correlation between methods
cat("\nCorrelation between LogNorm and scran size factors:\n")
cat("Overall Pearson r =",
    round(cor(sf_df$LogNorm_SF, sf_df$Scran_SF), 3), "\n")

# Do they agree differently for different cell types?
sf_cor_by_type <- sf_df %>%
  group_by(CellType) %>%
  summarise(r = round(cor(LogNorm_SF, Scran_SF), 3), .groups = "drop")
print(sf_cor_by_type)

# Scatter plot of size factors
p_sf <- sf_df %>%
  ggplot(aes(x = LogNorm_SF, y = Scran_SF, color = CellType)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("Neuron"    = "#4C72B0",
                                 "Astrocyte" = "#DD8452",
                                 "Microglia" = "#55A868")) +
  labs(title = "Size Factor Comparison: LogNormalize vs scran",
       subtitle = "Dashed line = perfect agreement",
       x = "LogNormalize Size Factor (= library size / 10K)",
       y = "scran Size Factor") +
  theme_bw(base_size = 13)

print(p_sf)

# ---- C2: UMAP comparison across all three methods ----

# Helper function to create UMAP plot
plot_umap <- function(seurat_obj_in, title, color_by = "true_type") {
  DimPlot(seurat_obj_in,
          reduction  = "umap",
          group.by   = color_by,
          cols       = c("Neuron"    = "#4C72B0",
                         "Astrocyte" = "#DD8452",
                         "Microglia" = "#55A868"),
          pt.size    = 0.8,
          label      = TRUE,
          label.size = 4) +
    labs(title = title) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

umap_log   <- plot_umap(seurat_log,   "LogNormalize")
umap_scran <- plot_umap(seurat_scran, "scran")
umap_sct   <- plot_umap(seurat_sct,   "SCTransform v2")

# Side-by-side comparison
combined_umap <- umap_log | umap_scran | umap_sct
print(combined_umap)

# ---- C3: Highly Variable Gene (HVG) overlap ----
# Which genes are selected as HVGs by each method?
hvg_log   <- VariableFeatures(seurat_log)
hvg_scran <- VariableFeatures(seurat_scran)
hvg_sct   <- VariableFeatures(seurat_sct)

# Are the true marker genes being captured?
# (Recall: markers are genes 1-50 for Neuron, 51-100 Astrocyte, 101-150 Microglia)
true_markers <- paste0("Gene", 1:150)

hvg_overlap <- data.frame(
  Method     = c("LogNormalize", "scran", "SCTransform"),
  Total_HVG  = c(length(hvg_log), length(hvg_scran), length(hvg_sct)),
  Markers_Recovered = c(
    sum(true_markers %in% hvg_log),
    sum(true_markers %in% hvg_scran),
    sum(true_markers %in% hvg_sct)
  )
)
hvg_overlap$Recall_Pct <- round(
  100 * hvg_overlap$Markers_Recovered / length(true_markers), 1)

cat("\n=== Highly Variable Gene Selection (True Markers Recovered) ===\n")
print(hvg_overlap)

# ---- C4: Marker gene expression comparison ----
# Plot expression of a Microglia-specific marker gene
# across normalization methods — does it look different?

microglia_marker_gene <- "Gene101"  # known microglia marker in our simulation

# Extract normalized expression from each method
expr_log   <- GetAssayData(seurat_log,   assay = "RNA", layer = "data")[microglia_marker_gene, ]
expr_scran <- GetAssayData(seurat_scran, assay = "RNA", layer = "data")[microglia_marker_gene, ]
expr_sct   <- GetAssayData(seurat_sct,   assay = "SCT", layer = "scale.data")[microglia_marker_gene, ]

expr_compare <- data.frame(
  CellType   = seurat_obj$true_type,
  LogNorm    = expr_log[colnames(seurat_obj)],
  Scran      = expr_scran[colnames(seurat_obj)],
  SCTransform = expr_sct[colnames(seurat_obj)]
)

# Reshape for plotting
expr_long <- tidyr::pivot_longer(expr_compare,
                                  cols      = c("LogNorm", "Scran", "SCTransform"),
                                  names_to  = "Method",
                                  values_to = "Expression")

p_marker <- ggplot(expr_long, aes(x = CellType, y = Expression, fill = CellType)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  facet_wrap(~Method, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Neuron"    = "#4C72B0",
                                "Astrocyte" = "#DD8452",
                                "Microglia" = "#55A868")) +
  labs(title = paste("Expression of", microglia_marker_gene, "(Microglia Marker)"),
       subtitle = "By normalization method — does each method correctly identify microglia-specific expression?",
       x = "Cell Type", y = "Normalized Expression") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "gray90"))

print(p_marker)
```

---

### Part D: Interpretation Questions

Answer the following questions in your lab notebook. 

> **Question 1 — Library Size Distribution**
>
> Look at the violin plot from Part B (`p_lib`). Microglia have substantially higher library sizes than Neurons and Astrocytes.
>
> **(a)** If you used `LogNormalize` (which divides all cells by their library size), what would happen to the expression of microglia-specific marker genes relative to how they appear in the raw counts? Would microglia still cluster separately from other cell types?
>
> *Possible interpretation A:* Library size normalization would correctly remove the technical depth differences within each cell type, and since microglia markers are biologically up-regulated (not just due to higher total counts), microglia should still cluster separately.
>
> *Possible interpretation B:* Because all cells are normalized to the same scale, the *relative* increase of microglia markers is diluted — some of the "signal" driving microglia separation in PCA/UMAP may be due to their genuine higher total RNA, which `LogNormalize` removes. This could cause microglia to appear *less* distinct.
>
> *Possible interpretation C:* `LogNormalize` over-corrects for microglia library size. After normalization, microglia housekeeping genes appear lower than in Neurons/Astrocytes (because the divisor is larger), potentially making microglia appear *more* similar to the other types than they truly are.
>
> *Which interpretation(s) are supported by your UMAP results? What additional analysis could you do to distinguish between interpretations B and C?*

---

> **Question 2 — Size Factor Comparison**
>
> Examine the scatter plot comparing LogNormalize and scran size factors (`p_sf`).
>
> **(a)** If the two methods agreed perfectly, all points would fall on the dashed diagonal line. Which cell type shows the most deviation? What does this tell you about how the two methods treat this cell type differently?
>
> **(b)** scran size factors should be roughly proportional to library size within a cell type, but can differ *between* cell types if pre-clustering properly identifies them. Explain why this is desirable behavior.
>
> **(c)** In what biological scenario would `scran` give **worse** results than `LogNormalize`? (Hint: think about the pre-clustering step.)

---

> **Question 3 — HVG Recall**
>
> The `hvg_overlap` table shows how many of the 150 true marker genes were recovered as highly variable genes by each method.
>
> **(a)** Which method recovers the most marker genes? Does this surprise you?
>
> **(b)** Would recovering **more** marker genes always be better? Consider a scenario where you are trying to identify a rare cell type (1% of cells) — what properties of normalization would you prioritize?
>
> **(c)** SCTransform's Pearson residuals are placed in the `scale.data` slot, not the `data` slot. What does this mean for downstream steps like differential expression testing? (Hint: look at the Seurat documentation for `FindMarkers` and which slot it uses by default for SCT.)

---

> **Question 4 — The "Right" Answer**
>
> In this simulation, we know the **ground truth** (true cell type labels, true marker genes). In a real experiment, you would not have this.
>
> **(a)** What criteria would you use in a real experiment to evaluate whether your normalization was successful?
>
> **(b)** A colleague argues: "All three methods give similar UMAPs, so it doesn't matter which one you use." Do you agree? Under what conditions would the choice of normalization method matter most?

---

### Part E: Code Modification Challenge

Choose **one** of the following modifications to implement. Document your changes with comments, run the modified analysis, and describe how (and why) the results differ from the original.

**Option 1 — Regress out confounders with SCTransform**

The original `SCTransform` call does not regress any additional covariates. In real datasets, you often want to regress out the fraction of reads mapping to mitochondrial genes (`percent.mt`), which is a proxy for cell stress/damage.

Modify the `SCTransform` call to:
1. Compute a simulated `percent.mt` variable (e.g., `runif(ncol(seurat_obj), 0, 0.3)` for demonstration)
2. Add it to the Seurat object metadata
3. Include `vars.to.regress = "percent.mt"` in `SCTransform`
4. Re-run PCA and UMAP, and compare to the original SCTransform UMAP

*Does regressing out a random variable change the clustering? What would you conclude from this?*

---

**Option 2 — Explore the effect of scale factor in LogNormalize**

The default scale factor is 10,000. Modify the analysis to run `LogNormalize` with scale factors of 100, 1,000, 10,000, and 100,000. 

1. Create a UMAP for each scale factor
2. Compare the resulting UMAPs and PCA variance explained
3. Compare the HVG lists

*Hint:* The log1p transformation means the scale factor affects the magnitude but not the rank of normalized values. What does this imply about its practical effect on clustering?

---

**Option 3 — Introduce a batch effect and test correction**

Simulate a simple batch effect by adding a constant offset to half the cells' counts:

```r
# After creating counts_matrix_dropout but before making the Seurat object:
batch <- c(rep("Batch1", 250), rep("Batch2", 250))

# Add a systematic batch offset to Batch 2 cells
batch2_cells <- which(batch == "Batch2")
counts_batch <- counts_matrix_dropout
counts_batch[, batch2_cells] <- counts_batch[, batch2_cells] + 
  matrix(rpois(nrow(counts_batch) * length(batch2_cells), lambda = 3),
         nrow = nrow(counts_batch))

# Proceed with the modified count matrix...
```

After applying all three normalization methods to this batch-affected data:
1. Color UMAPs by `batch` instead of `true_type` — does the batch drive clustering?
2. Apply `Seurat::IntegrateData` or `harmony` after normalization to correct the batch
3. Compare marker gene recall before and after batch correction

*Does normalization alone remove batch effects? When is dedicated batch correction necessary?*

---

### Key Takeaways

1. **No single normalization method is universally correct.** The best choice depends on your biological question, cell type composition, and data characteristics.

2. **Normalization assumptions matter.** `LogNormalize` assumes equal total RNA per cell. `scran` is more flexible but requires pre-clustering. `SCTransform` makes no assumption about equal totals but produces residuals rather than counts.

3. **Always visualize your data before and after normalization.** Library size distributions, mean-variance relationships, and technical covariate distributions should all be examined.

4. **Normalization is not a panacea.** It cannot remove batch effects, doublets, or ambient RNA contamination — these require dedicated QC steps.

5. **Document your normalization choices.** A reproducible analysis must specify the exact method, parameters, package versions, and random seed used. Use `sessionInfo()` at the end of every analysis script.

---

*Module prepared for BIOINF 504: Rigor & Reproducibility in Bioinformatics. Please report any errors or suggestions to the course instructor.*

*Last updated: 2026 | R packages: Seurat ≥ 5.0, scran ≥ 1.28, scater ≥ 1.28*

---

## Key References

- Mortazavi A et al. (2008). Mapping and quantifying mammalian transcriptomes by RNA-Seq. *Nature Methods*, 5(7):621–628. *(RPKM)*
- Wagner GP, Kin K, Lynch VJ (2012). Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. *Theory in Biosciences*, 131(4):281–285. *(RPKM vs TPM)*
- Li B & Dewey CN (2011). RSEM: accurate transcript quantification from RNA-seq data with or without a reference genome. *BMC Bioinformatics*, 12:323. *(TPM)*
- Zhao S, Ye Z, Stanton R (2020). Misuse of RPKM or TPM normalization when comparing across samples and sequencing protocols. *RNA*, 26(8):903–909. *(TPM limitations)*
- Robinson MD & Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*, 11(3):R25. *(TMM)*
- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12):550. *(DESeq2)*
- Lun ATL, Bach K, Marioni JC (2016). Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. *Genome Biology*, 17:75. *(scran)*
- Hafemeister C & Satija R (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology*, 20:296. *(SCTransform)*
- Choudhary S & Satija R (2022). Comparison and evaluation of statistical error models for scRNA-seq. *Genome Biology*, 23:27. *(SCTransform v2)*
- Huber W et al. (2002). Variance stabilization applied to microarray data calibration. *Bioinformatics*, 18(Suppl 1):S96–S104. *(VSN)*
- Dieterle F et al. (2006). Probabilistic quotient normalization as robust method to account for dilution of complex biological mixtures. *Analytical Chemistry*, 78(13):4281–4290. *(PQN)*
