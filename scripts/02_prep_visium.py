"""
02_prep_visium.py — Load the 6 Wu et al. Visium sections, attach pathology ground
truth + spot coords, and export a clean R handoff (the R<->Python boundary).

GOTCHA handled: Wu's *.gz files are actually PLAINTEXT (misleading extension), so we
read the MatrixMarket via an open binary stream (bypassing scipy's extension-based
gunzip) and re-export as REAL gzip for R.

Outputs:
  data/processed/visium/<section>.h5ad                  # spots x genes, for Part B / Python validation
  data/processed/rctd_input/<section>/{counts.mtx.gz,genes.tsv,barcodes.tsv,spot_meta.csv}
  data/processed/rctd_input/reference/{counts.mtx.gz,genes.tsv,cell_meta.csv}   # shared scRNA reference
"""
import os, scipy.io, scipy.sparse as sp, numpy as np, pandas as pd, anndata as ad, scanpy as sc, gzip

RAW = "/home/ubuntu/spatial-tme-pipeline/data/raw"
PROC = "/home/ubuntu/spatial-tme-pipeline/data/processed"
REF = f"{PROC}/reference_annotated.h5ad"
SECS = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]


def read_plaintext_mtx(path):          # file is plaintext despite .gz
    with open(path, "rb") as fh:
        return scipy.io.mmread(fh).tocsr()


def load_lines(path):
    with open(path) as fh:
        return [l.strip().split("\t")[0] for l in fh if l.strip()]


def write_mtx_gz(mat, path):           # write genes x cells as REAL gzip MatrixMarket
    with gzip.open(path, "wb") as fh:
        scipy.io.mmwrite(fh, sp.csr_matrix(mat))


def export_reference():
    a = sc.read_h5ad(REF)              # cells x genes, X = counts
    d = f"{PROC}/rctd_input/reference"; os.makedirs(d, exist_ok=True)
    write_mtx_gz(a.X.T, f"{d}/counts.mtx.gz")                       # genes x cells
    pd.Series(a.var_names).to_csv(f"{d}/genes.tsv", index=False, header=False)
    a.obs[["cell_type", "compartment", "sample", "subtype"]].to_csv(f"{d}/cell_meta.csv")
    # marker genes per cell type for SPOTlight (mgs: gene/cluster/weight)
    b = a.copy(); sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    sc.tl.rank_genes_groups(b, "cell_type", method="wilcoxon")
    rows = []
    for ct in b.obs["cell_type"].cat.categories:
        df = sc.get.rank_genes_groups_df(b, group=ct).head(100)
        df = df[df["logfoldchanges"] > 0]
        rows.append(pd.DataFrame({"gene": df["names"], "cluster": ct,
                                  "weight": df["logfoldchanges"].clip(lower=0)}))
    pd.concat(rows).to_csv(f"{d}/mgs.csv", index=False)
    print(f"reference exported: {a.shape[0]} cells x {a.shape[1]} genes + mgs -> {d}")


def load_section(sec):
    mdir = f"{RAW}/filtered_count_matrices/{sec}_filtered_count_matrix"
    X = read_plaintext_mtx(f"{mdir}/matrix.mtx.gz")                 # genes x spots
    genes = load_lines(f"{mdir}/features.tsv.gz")
    barcodes = load_lines(f"{mdir}/barcodes.tsv.gz")
    a = ad.AnnData(sp.csr_matrix(X.T), dtype=np.float32)           # spots x genes
    a.var_names = pd.Index(genes).astype(str); a.var_names_make_unique()
    a.obs_names = barcodes

    md = pd.read_csv(f"{RAW}/metadata/{sec}_metadata.csv", index_col=0)
    pos_p = f"{RAW}/spatial/{sec}_spatial/tissue_positions_list.csv"
    pos = pd.read_csv(pos_p, header=None,
                      names=["barcode", "in_tissue", "array_row", "array_col", "px_row", "px_col"]).set_index("barcode")

    keep = [b for b in a.obs_names if b in md.index]               # annotated in-tissue spots
    a = a[keep].copy()
    a.obs["pathology"] = md.loc[keep, "Classification"].values
    a.obs["subtype"] = md.loc[keep, "subtype"].values
    a.obs["section"] = sec
    xy = pos.reindex(keep)[["px_col", "px_row"]].values.astype(float)
    a.obsm["spatial"] = xy
    return a


def export_section(a, sec):
    d = f"{PROC}/rctd_input/{sec}"; os.makedirs(d, exist_ok=True)
    write_mtx_gz(a.X.T, f"{d}/counts.mtx.gz")                      # genes x spots
    pd.Series(a.var_names).to_csv(f"{d}/genes.tsv", index=False, header=False)
    pd.Series(a.obs_names).to_csv(f"{d}/barcodes.tsv", index=False, header=False)
    meta = a.obs[["pathology", "subtype", "section"]].copy()
    meta["x"] = a.obsm["spatial"][:, 0]; meta["y"] = a.obsm["spatial"][:, 1]
    meta.to_csv(f"{d}/spot_meta.csv")


def main():
    os.makedirs(f"{PROC}/visium", exist_ok=True)
    export_reference()
    for sec in SECS:
        a = load_section(sec)
        a.write(f"{PROC}/visium/{sec}.h5ad")
        export_section(a, sec)
        xmax = float(a.X.max())
        print(f"{sec}: {a.n_obs} spots x {a.n_vars} genes | counts_max={xmax:.0f} | "
              f"pathology cats={a.obs['pathology'].nunique()}")


if __name__ == "__main__":
    main()
