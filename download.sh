#!/usr/bin/env bash
# download.sh — fetch and organize all input data for the spatial-tme-pipeline.
#
# PROVENANCE
#   Spatial (Visium):  Wu et al. 2021, Nature Genetics — "A single-cell and spatially
#                      resolved atlas of human breast cancers." Zenodo record 4739739
#                      (DOI 10.5281/zenodo.4739739). 6 primary breast cancer sections
#                      with per-spot pathologist annotation (the validation ground truth).
#   scRNA reference:   GSE161529 (Pal et al.) — 2 ER+ + 2 TNBC tumors, 31,265 QC'd cells.
#                      Produced by the companion repo `scrna-tme-pipeline` (bulk->scRNA->
#                      spatial through-line). See the REFERENCE section at the bottom.
#
# NOTE: Wu's *.gz files are plaintext despite the extension — scripts/02_prep_visium.py
#       handles that; nothing to do here.
set -euo pipefail
RAW="data/raw"; mkdir -p "$RAW"

# ---------------------------------------------------------------------------
# 1. Wu et al. Visium data (Zenodo 4739739)
# ---------------------------------------------------------------------------
ZEN="https://zenodo.org/records/4739739/files"
for f in filtered_count_matrices spatial metadata; do   # (skip raw_count_matrices + images.pdf: not needed)
  if [ -d "$RAW/$f" ]; then echo "[skip] $RAW/$f exists"; continue; fi
  echo "[get ] $f.tar.gz"
  wget -q -O "$RAW/$f.tar.gz" "$ZEN/$f.tar.gz?download=1"
  tar -xzf "$RAW/$f.tar.gz" -C "$RAW"
  rm -f "$RAW/$f.tar.gz"
done
echo "Visium ready: 6 sections in $RAW/{filtered_count_matrices,spatial,metadata}/"

# ---------------------------------------------------------------------------
# 2. scRNA reference (GSE161529)
# ---------------------------------------------------------------------------
# The deconvolution reference (data/processed/reference_annotated.h5ad) is built by
# scripts/01_build_reference.py from the annotated single-cell object of the companion
# project. To regenerate that input:
#
#   git clone https://github.com/pipelinepilot1999/scrna-tme-pipeline.git ~/scrna-tme-pipeline
#   cd ~/scrna-tme-pipeline && bash download.sh          # pulls the 4 GSM samples from GEO
#   # run its QC notebook -> data/processed/01_qc_filtered.h5ad  (31,265 cells, raw counts)
#
# Then, back in this repo:
#   python scripts/01_build_reference.py                 # marker-based 7-type annotation
#
# GEO GSE161529 samples used (2 ER+, 2 TNBC):
#   GSM4909296 ER-MH0001 | GSM4909301 ER-MH0042 | GSM4909281 TN-MH0126 | GSM4909282 TN-MH0135
echo "scRNA reference: see section 2 above (companion repo scrna-tme-pipeline / GSE161529)."
