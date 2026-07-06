# Spatial Deconvolution of the Breast Cancer Tumor Microenvironment

**Reference-based deconvolution of breast cancer Visium sections, validated against
pathologist annotations.** The capstone of a three-part breast-cancer through-line.


![Deconvolved composition of 6 breast cancer Visium sections](figures/hero_composition.png)
*Each Visium spot colored by its RCTD-deconvolved composition — R = malignant, G = immune,
B = stromal. The spatial "where" that bulk and dissociated single-cell cannot resolve.*

## The through-line: bulk → single-cell → spatial

| Stage | Question it answers | What it cannot do |
|-------|--------------------|-------------------|
| **Bulk RNA-seq** | Is there an ER+/TNBC signature? (GATA3/ESR1 vs KRT14/FOXC1) | Averages all cells — can't say *which* cells |
| **scRNA-seq** (GSE161529) | *Which cells* carry it? Are they malignant? | Localizes to malignant epithelium, proves malignancy cell-by-cell via CNV — but tissue is dissociated, spatial context destroyed |
| **Spatial (this repo)** | *Where* do malignant / immune / stromal populations sit relative to each other? | Puts the compartments back into tissue context |

Bulk found the signature by averaging; scRNA localized it to malignant epithelium and
proved malignancy by CNV; **spatial restores the "where"** — the resolution neither
bulk nor dissociated single-cell can provide. Part B (niche analysis with squidpy)
builds on these proportions to characterize the tumor–immune boundary by subtype.

## What this does (Part A)

1. Rebuilds an annotated **scRNA reference** (GSE161529 / Pal et al.; 31,265 cells, 7
   cell types) from raw counts.
2. **Deconvolves** each Wu et al. Visium section into per-spot cell-type **proportion
   vectors** using two independent methods.
3. **Validates** the proportions against the pathologist's per-spot annotations —
   reporting concordance *and* where it fails.

## Method choice (and why)

**Primary: RCTD (`spacexr`, Cable et al. 2022).** Purpose-built for spatial
deconvolution, CPU-native, returns per-spot proportion weights, and includes an
explicit platform-effect normalization for the scRNA→Visium capture difference. Verified
running on this CPU-only r5 instance (`max_cores=4`; the `max_cores=1` `fitBulk` hang is
avoided). Runs in "full" mode → a proportion vector per spot.

**Cross-check: SPOTlight (NMF).** An independent method with a different statistical
basis (seeded NMF vs RCTD's Poisson regression). "Two independent methods agree" is the
rigor signal; disagreement is reported, not smoothed over.

*Not used:* cell2location (field-standard but Bayesian/GPU-oriented — hours on CPU here).

### The R ↔ Python boundary (a deliberate design choice)
Deconvolution (RCTD/SPOTlight) is **R**; the reference build, prep, validation, and
Part B (squidpy) are **Python**. The handoff is explicit and file-based: Python exports
counts + coords + barcodes (MatrixMarket + CSV); R writes per-spot proportion CSVs;
Python reads them back. No shared-object gymnastics, fully reproducible across the
language boundary.

## Reference annotation — an honest note
The original scRNA repo's clustering notebook had lost its Harmony step, and its cell
types were keyed to *hardcoded leiden cluster numbers* (unstable across re-runs). Rather
than risk a silent mis-label, the reference here reproduces the same pipeline
(QC→HVG→PCA→Harmony→leiden 0.5) and the same 7 cell types, but assigns each cluster by
**canonical marker enrichment** (per-cluster mean marker expression, z-scored per gene
across clusters — a self-validating, renumbering-robust method). The per-cluster evidence
is saved to `data/processed/reference_cluster_evidence.csv`.

Cell types: `ER+ tumor`, `TNBC tumor`, `T cell`, `Myeloid`, `Plasma cell`, `Fibroblast`,
`Endothelial` (grouped into malignant / immune / stromal compartments for validation).

## Data
- **Reference:** GSE161529 (Pal et al.) — 2 ER+ + 2 TNBC tumors, 31,265 QC'd cells.
- **Spatial:** Wu et al. 2021, *Nat Genet* — Visium, **Zenodo 10.5281/zenodo.4739739**
  (note: 4739739, not the 4739749 sometimes cited). 6 sections (2 ER+, 4 TNBC),
  1.1k–4.9k spots each, with **per-spot pathologist `Classification`** = the independent
  ground truth.

## Validation
Deconvolved all 6 sections (15,608 pathologist-annotated spots) into 7-cell-type
proportion vectors with RCTD, cross-checked with SPOTlight.

**Malignant signal tracks the pathology.** Spots the pathologist called cancer-containing
carry a median malignant proportion of **0.62**, vs **0.37** for non-cancer spots
(Mann–Whitney p≈0; AUROC **0.75** for malignant-proportion as a cancer classifier). RCTD
resolves the *subtype* correctly per section: TNBC sections are TNBC-tumor-dominant
(1142243F: TNBC 0.55 / ER+ 0.02), ER+ sections are ER+-tumor-dominant (CID4535: ER+ 0.60 /
TNBC 0.01). Pure "Lymphocytes" spots carry a median immune proportion of **0.58**.

**Two independent methods agree.** RCTD vs SPOTlight malignant proportion correlates at
**r = 0.87–0.97 per section** (80% of spots share the same dominant compartment). SPOTlight
corroborates *where* the malignant compartment is, though it does not resolve the ER+/TNBC
subtype that RCTD recovers, and is noisier — a known SPOTlight accuracy gap.
(`figures/rctd_vs_spotlight.png`)

**Where it fails — reported, not hidden:**
- *Stroma is underestimated.* Pure-"Stroma" spots carry only a median 0.24 stromal
  proportion, and stroma-labeled spots are mostly predicted malignant/immune. Two causes:
  the reference is only ~5% stromal cells (1,448 fibroblast + 128 endothelial of 31,265),
  and breast-tumor stroma is genuinely tumor-infiltrated at spot resolution. This drives the
  modest overall dominant-compartment accuracy (0.40 on pure-label spots); malignant (87%)
  and immune (78%) recovery are good, stroma (17%) is the weak point.
- *No normal-epithelial reference type.* Pathologist "normal gland/duct" spots get a median
  malignant proportion of 0.56 — normal epithelium is forced into the malignant profile
  because the reference has no normal-epithelial class. A structural limitation, quantified.

Figures: `hero_composition.png` (tissue colored by composition), `composition_vs_pathology.png`,
`malignant_by_pathology.png`, `rctd_vs_spotlight.png`. Full log: `data/processed/validation_summary.txt`.


![Malignant proportion by pathology label](figures/malignant_by_pathology.png)
![RCTD composition vs pathologist annotation](figures/composition_vs_pathology.png)
![RCTD vs SPOTlight cross-check](figures/rctd_vs_spotlight.png)

## Pipeline
```
scripts/01_build_reference.py   # scRNA reference (raw counts + 7 marker-based cell types)
scripts/02_prep_visium.py       # load 6 Visium sections + pathology + coords; export R handoff
scripts/03_run_rctd.R           # RCTD deconvolution (primary)  -> per-spot proportions
scripts/04_run_spotlight.R      # SPOTlight cross-check          -> per-spot proportions
scripts/05_validate.py          # concordance vs pathology + honest failure probes
scripts/06_figures.py           # hero composition figure
```

## Environments
- `spatial-r` — R 4.5.3, `bioconductor-spacexr 1.2.0`, `bioconductor-spotlight 1.14.0`
- `scrna` / `spatial` — Python (scanpy, anndata; squidpy for Part B)

## Limitations (read this)
- Deconvolution returns **proportions, not ground-truth cell identity** — a spot's
  malignant proportion is a mixture estimate, not a cell count.
- The reference has **no normal-epithelial type**, so pathologist "normal gland/duct"
  spots are forced toward the malignant profile — a known, quantified failure mode.
- **Malignant cells are patient-specific**; a reference from 4 patients may transfer
  imperfectly to other patients' tumors (a general limitation of scRNA→spatial transfer).
- RCTD's platform normalization mitigates but does not perfectly remove the scRNA↔Visium
  capture difference.
- Single cohort, 6 sections — scope is deliberately bounded.

## References
- Cable et al. *Robust decomposition of cell type mixtures in spatial transcriptomics.* Nat Biotechnol 2022. (RCTD)
- Wu et al. *A single-cell and spatially resolved atlas of human breast cancers.* Nat Genet 2021. (Visium data)
- Pal et al. GSE161529. (scRNA reference)
