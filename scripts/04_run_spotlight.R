# 04_run_spotlight.R — SPOTlight (NMF) deconvolution, INDEPENDENT CROSS-CHECK.
# Second method with a different statistical basis (seeded NMF vs RCTD's Poisson
# regression). Agreement between the two is the honest-rigor signal. Same handoff.
suppressMessages({library(SPOTlight); library(Matrix)})

PROC <- "/home/ubuntu/spatial-tme-pipeline/data/processed"
OUTD <- file.path(PROC, "spotlight_output"); dir.create(OUTD, showWarnings = FALSE)
SECS <- c("1142243F","1160920F","CID4290","CID4465","CID44971","CID4535")

read_counts <- function(dir) {
  m <- as(readMM(gzfile(file.path(dir, "counts.mtx.gz"))), "CsparseMatrix")
  rownames(m) <- readLines(file.path(dir, "genes.tsv")); m
}

refdir <- file.path(PROC, "rctd_input", "reference")
refm   <- read_counts(refdir)
cellmeta <- read.csv(file.path(refdir, "cell_meta.csv"), row.names = 1, check.names = FALSE)
colnames(refm) <- rownames(cellmeta)
groups <- as.character(cellmeta$cell_type)
mgs <- read.csv(file.path(refdir, "mgs.csv"), stringsAsFactors = FALSE)
cat("reference:", nrow(refm), "x", ncol(refm), "| mgs rows:", nrow(mgs), "\n")

for (sec in SECS) {
  cat("\n==== SPOTlight", sec, "====\n")
  secdir <- file.path(PROC, "rctd_input", sec)
  spm <- read_counts(secdir); colnames(spm) <- readLines(file.path(secdir, "barcodes.tsv"))
  res <- SPOTlight(x = refm, y = spm, groups = groups, mgs = mgs,
                   gene_id = "gene", group_id = "cluster", weight_id = "weight",
                   slot_sc = "counts", slot_sp = "counts", verbose = FALSE)
  props <- res$mat                                   # spots x cell_types
  outp <- file.path(OUTD, paste0(sec, "_spotlight_proportions.csv"))
  write.csv(round(props, 5), outp)
  cat("wrote", outp, "| dim", nrow(props), "x", ncol(props), "\n")
}
cat("\nALL_SPOTLIGHT_DONE\n")
