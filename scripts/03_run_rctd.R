# 03_run_rctd.R — RCTD (spacexr) deconvolution, PRIMARY method.
# Reads the CSV/MTX handoff from Python (02_prep_visium.py), runs RCTD full mode
# (per-spot proportion vectors) on each Wu et al. section, writes proportions CSV.
#
# spacexr 1.2.0 (Bioconductor) S4 API: createRctd(spatial_spe, reference_se) -> runRctd().
# Full mode output = assay(res,"weights") as cell_types x spots; we transpose to spots x types.
suppressMessages({library(spacexr); library(SpatialExperiment); library(SummarizedExperiment); library(Matrix)})

PROC <- "/home/ubuntu/spatial-tme-pipeline/data/processed"
OUTD <- file.path(PROC, "rctd_output"); dir.create(OUTD, showWarnings = FALSE)
SECS <- c("1142243F","1160920F","CID4290","CID4465","CID44971","CID4535")
NCORES <- 4

read_counts <- function(dir) {           # genes x cells/spots sparse matrix
  m <- as(readMM(gzfile(file.path(dir, "counts.mtx.gz"))), "CsparseMatrix")
  rownames(m) <- readLines(file.path(dir, "genes.tsv"))
  m
}

# --- reference (shared) ---
refdir <- file.path(PROC, "rctd_input", "reference")
refm   <- read_counts(refdir)
cellmeta <- read.csv(file.path(refdir, "cell_meta.csv"), row.names = 1, check.names = FALSE)
colnames(refm) <- rownames(cellmeta)
cellmeta$cell_type <- factor(cellmeta$cell_type)
reference_se <- SummarizedExperiment(assays = list(counts = refm),
                                     colData = DataFrame(cell_type = cellmeta$cell_type))
cat("reference:", nrow(refm), "genes x", ncol(refm), "cells |",
    nlevels(cellmeta$cell_type), "types\n")

for (sec in SECS) {
  cat("\n==== RCTD", sec, "====\n")
  secdir <- file.path(PROC, "rctd_input", sec)
  spm  <- read_counts(secdir)
  bc   <- readLines(file.path(secdir, "barcodes.tsv"))
  meta <- read.csv(file.path(secdir, "spot_meta.csv"), row.names = 1, check.names = FALSE)
  colnames(spm) <- bc
  coords <- as.matrix(meta[bc, c("x","y")])
  spatial_spe <- SpatialExperiment(assay = list(counts = spm), spatialCoords = coords)

  rctd <- createRctd(spatial_spe, reference_se, cell_type_col = "cell_type")
  res  <- runRctd(rctd, rctd_mode = "full", max_cores = NCORES)
  w <- as.matrix(assay(res, "weights"))          # cell_types x spots
  props <- t(w)                                   # spots x cell_types
  props <- props / rowSums(props)                 # normalize to proportions
  outp <- file.path(OUTD, paste0(sec, "_rctd_proportions.csv"))
  write.csv(round(props, 5), outp)
  cat("wrote", outp, "| dim", nrow(props), "x", ncol(props), "\n")
}
cat("\nALL_RCTD_DONE\n")
