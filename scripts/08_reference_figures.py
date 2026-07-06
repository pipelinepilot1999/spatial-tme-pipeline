"""
08_reference_figures.py — Visual evidence for the scRNA reference annotation.
UMAP of the 31k cells colored by the 7 assigned cell types, and a canonical-marker
dotplot showing the basis for each label (the annotation is marker-based, so this is
the audit figure). Recomputes the Harmony embedding used for clustering.
"""
import warnings, scanpy as sc, numpy as np, harmonypy, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
sc.settings.figdir = "/home/ubuntu/spatial-tme-pipeline/figures"
REF = "/home/ubuntu/spatial-tme-pipeline/data/processed/reference_annotated.h5ad"
MARK = {"ER+ tumor":["ESR1","GATA3","FOXA1","TFF3"], "TNBC tumor":["KRT14","KRT17","KRT5","FOXC1"],
        "T cell":["CD3D","CD3E","CD8A"], "Myeloid":["CD68","LYZ","TYROBP"],
        "Plasma cell":["MZB1","JCHAIN","IGHG1"], "Fibroblast":["COL1A1","DCN","LUM"],
        "Endothelial":["PECAM1","VWF","CLDN5"]}
ORD = list(MARK)

a = sc.read_h5ad(REF)                       # X = counts
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
a.raw = a
hv = a.copy(); sc.pp.highly_variable_genes(hv, n_top_genes=2000, flavor="seurat")
hv = hv[:, hv.var.highly_variable].copy(); sc.pp.scale(hv, max_value=10)
sc.pp.pca(hv, n_comps=50, random_state=0)
ho = harmonypy.run_harmony(hv.obsm["X_pca"], hv.obs, ["sample"], random_state=0)
Z = ho.Z_corr; Z = Z.detach().cpu().numpy() if hasattr(Z,"detach") else np.asarray(Z)
if Z.shape[0] != hv.n_obs: Z = Z.T
a.obsm["X_pca_harmony"] = Z
sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony", random_state=0)
sc.tl.umap(a, random_state=0)
a.obs["cell_type"] = a.obs["cell_type"].astype("category")

sc.pl.umap(a, color="cell_type", title="scRNA reference — 7 cell types (GSE161529, 31,265 cells)",
           legend_fontsize=9, size=8, show=False, save=None)
plt.tight_layout(); plt.savefig(f"{sc.settings.figdir}/reference_umap_celltypes.png", dpi=140, bbox_inches="tight"); plt.close()
sc.pl.umap(a, color="subtype", title="scRNA reference — sample subtype", size=8, show=False)
plt.tight_layout(); plt.savefig(f"{sc.settings.figdir}/reference_umap_subtype.png", dpi=140, bbox_inches="tight"); plt.close()

sc.pl.dotplot(a, MARK, groupby="cell_type", categories_order=ORD, standard_scale="var",
              title="Canonical markers by assigned cell type (annotation basis)", show=False)
plt.savefig(f"{sc.settings.figdir}/reference_marker_dotplot.png", dpi=140, bbox_inches="tight"); plt.close()
print("wrote reference_umap_celltypes.png, reference_umap_subtype.png, reference_marker_dotplot.png")
print(a.obs["cell_type"].value_counts())
