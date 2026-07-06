"""
01_build_reference.py — Regenerate the annotated scRNA reference for deconvolution.

Source: ~/scrna-tme-pipeline/data/processed/01_qc_filtered.h5ad
        (31,265 QC'd cells, RAW COUNTS, all genes; GSE161529 / Pal et al.).

Reproduces Sai's pipeline (QC->normalize->HVG->PCA->Harmony->leiden 0.5), then labels
each cluster by CANONICAL MARKER ENRICHMENT: mean marker expression per cluster,
z-scored PER GENE across clusters (so cross-cell-type comparison is on one scale),
averaged per signature, argmax = cell type. Epithelial clusters split ER+ vs TNBC by
ER-program vs basal-program z. Self-validating + full audit trail logged.

Output: data/processed/reference_annotated.h5ad
  X=raw counts (all genes); obs['cell_type'] (7 types), obs['compartment']
  (malignant/immune/stromal), obs['sample','subtype'].
Also: data/processed/reference_cluster_evidence.csv
"""
import scanpy as sc, pandas as pd, numpy as np, harmonypy

SRC  = "/home/ubuntu/scrna-tme-pipeline/data/processed/01_qc_filtered.h5ad"
OUT  = "/home/ubuntu/spatial-tme-pipeline/data/processed/reference_annotated.h5ad"
EVID = "/home/ubuntu/spatial-tme-pipeline/data/processed/reference_cluster_evidence.csv"
sc.settings.verbosity = 1

# Canonical markers (Sai's sets, broadened with textbook markers for robustness).
SIG = {
    "Epithelial":  ["EPCAM","KRT8","KRT18","KRT19","ELF3","TACSTD2","KRT5","KRT14","KRT17"],
    "T cell":      ["CD3D","CD3E","CD8A","CD2","TRAC"],
    "Myeloid":     ["CD68","LYZ","CD14","TYROBP","FCER1G","AIF1"],
    "Plasma cell": ["MZB1","JCHAIN","IGHG1","DERL3","XBP1"],
    "Fibroblast":  ["COL1A1","COL1A2","DCN","LUM","PDGFRB"],
    "Endothelial": ["PECAM1","VWF","CLDN5","CD34","EGFL7"],
}
ER_PROG    = ["ESR1","GATA3","FOXA1","AR","TFF3"]
BASAL_PROG = ["KRT5","KRT14","KRT17","FOXC1","MIA"]
COMPART = {"ER+ tumor":"malignant","TNBC tumor":"malignant","T cell":"immune",
           "Myeloid":"immune","Plasma cell":"immune","Fibroblast":"stromal","Endothelial":"stromal"}


def cluster_marker_z(adata, genes_by_sig):
    """mean expr per cluster per marker (log-norm), z-scored per gene across clusters."""
    allg = sorted({g for gs in genes_by_sig.values() for g in gs if g in adata.raw.var_names})
    idx = [adata.raw.var_names.get_loc(g) for g in allg]
    X = adata.raw.X[:, idx]
    df = pd.DataFrame(X.toarray() if hasattr(X, "toarray") else np.asarray(X),
                      columns=allg, index=adata.obs_names)
    df["leiden"] = adata.obs["leiden"].values
    cmean = df.groupby("leiden", observed=True).mean()
    z = (cmean - cmean.mean(0)) / (cmean.std(0) + 1e-9)         # clusters x genes
    sigz = pd.DataFrame({s: z[[g for g in gs if g in z.columns]].mean(1)
                         for s, gs in genes_by_sig.items()})
    return sigz


def main():
    adata = sc.read_h5ad(SRC)
    print("loaded:", adata.shape, "| obs:", list(adata.obs.columns))
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    hvg = adata.copy()
    sc.pp.highly_variable_genes(hvg, n_top_genes=2000, flavor="seurat")
    hvg = hvg[:, hvg.var.highly_variable].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.pp.pca(hvg, n_comps=50, random_state=0)
    ho = harmonypy.run_harmony(hvg.obsm["X_pca"], hvg.obs, ["sample"], random_state=0)
    Z = ho.Z_corr
    if hasattr(Z, "detach"): Z = Z.detach().cpu().numpy()
    Z = np.asarray(Z)
    if Z.shape[0] != hvg.n_obs: Z = Z.T
    hvg.obsm["X_pca_harmony"] = Z
    sc.pp.neighbors(hvg, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony", random_state=0)
    sc.tl.leiden(hvg, resolution=0.5, flavor="leidenalg", random_state=0)
    adata.obs["leiden"] = hvg.obs["leiden"].values
    print("leiden clusters:", adata.obs["leiden"].nunique())

    sigz = cluster_marker_z(adata, {**SIG, "ER_prog": ER_PROG, "Basal_prog": BASAL_PROG})
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
    top = {cl: ", ".join(adata.uns["rank_genes_groups"]["names"][cl][:8]) for cl in sigz.index}
    subtypemix = adata.obs.groupby("leiden", observed=True)["subtype"].agg(
        lambda s: s.value_counts().idxmax())

    mapping = {}
    for cl in sigz.index:
        base = max(SIG, key=lambda k: sigz.loc[cl, k])
        if base == "Epithelial":
            base = "ER+ tumor" if sigz.loc[cl, "ER_prog"] >= sigz.loc[cl, "Basal_prog"] else "TNBC tumor"
        mapping[cl] = base
    adata.obs["cell_type"] = adata.obs["leiden"].map(mapping).astype("category")
    adata.obs["compartment"] = adata.obs["cell_type"].map(COMPART).astype("category")

    ev = sigz.round(2).copy()
    ev["n_cells"] = adata.obs["leiden"].value_counts()
    ev["dom_subtype"] = subtypemix
    ev["assigned"] = pd.Series(mapping)
    ev["top_genes"] = pd.Series(top)
    ev.to_csv(EVID)
    print("\n=== per-cluster marker z-scores + assignment ===")
    with pd.option_context("display.width", 240, "display.max_columns", 25):
        print(ev[["n_cells","dom_subtype","assigned","Epithelial","ER_prog","Basal_prog",
                  "T cell","Myeloid","Plasma cell","Fibroblast","Endothelial"]])
    for cl in ev.index:
        print(f"  cluster {cl:>2} [{mapping[cl]:<12}] top: {top[cl]}")

    adata.X = adata.layers["counts"]
    adata.obs = adata.obs[[c for c in adata.obs.columns]]
    del adata.raw
    print("\n=== cell_type breakdown ===\n", adata.obs["cell_type"].value_counts())
    print("\n=== compartment breakdown ===\n", adata.obs["compartment"].value_counts())
    xmax = float(adata.X.max())
    print(f"\nX is counts? max={xmax:.1f} integer-like={xmax.is_integer()}")
    adata.write(OUT)
    print("wrote", OUT, adata.shape)


if __name__ == "__main__":
    main()
