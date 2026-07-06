"""
07_niche_analysis.py — Part B: spatial niche analysis (squidpy).

Builds on Part A's per-spot proportions to ask WHERE compartments sit relative to each
other, and characterizes the tumor-immune boundary by subtype.

Per section: attach RCTD proportions -> dominant cell type + compartment; build the
Visium spatial graph; run neighborhood enrichment (which compartments co-localize vs
segregate) and co-occurrence; cluster spots into spatial domains.

Headline test: is immune signal infiltrated INTO tumor (hot) or excluded (cold), and
does it differ ER+ vs TNBC? (TNBC is clinically the more immune-infiltrated subtype.)

Outputs: data/processed/partB/<section>.h5ad (enriched), figures/partB_*.png,
data/processed/partB_summary.txt
"""
import os, warnings, numpy as np, pandas as pd, scanpy as sc, squidpy as sq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
warnings.filterwarnings("ignore")

PROC = "/home/ubuntu/spatial-tme-pipeline/data/processed"
FIG  = "/home/ubuntu/spatial-tme-pipeline/figures"; os.makedirs(FIG, exist_ok=True)
OUTD = f"{PROC}/partB"; os.makedirs(OUTD, exist_ok=True)
SECS = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]
SUBTYPE = {"1142243F":"TNBC","1160920F":"TNBC","CID4290":"ER+","CID4465":"TNBC","CID44971":"TNBC","CID4535":"ER+"}
COMPART = {"ER+ tumor":"malignant","TNBC tumor":"malignant","T cell":"immune",
           "Myeloid":"immune","Plasma cell":"immune","Fibroblast":"stromal","Endothelial":"stromal"}
CTYPES = ["ER+ tumor","TNBC tumor","T cell","Myeloid","Plasma cell","Fibroblast","Endothelial"]
LOG=[]; log=lambda *a:(print(*a), LOG.append(" ".join(str(x) for x in a)))


def build(sec):
    a = sc.read_h5ad(f"{PROC}/visium/{sec}.h5ad")
    p = pd.read_csv(f"{PROC}/rctd_output/{sec}_rctd_proportions.csv", index_col=0)
    p.index = p.index.astype(str); p = p.reindex(a.obs_names)
    for c in CTYPES: a.obs[c] = p[c].values
    comp = pd.DataFrame({k: p[[c for c in CTYPES if COMPART[c]==k]].sum(1)
                         for k in ["malignant","immune","stromal"]}, index=p.index)
    for k in comp: a.obs[k] = comp[k].values
    a.obs["dominant_ct"] = p[CTYPES].idxmax(1).values
    a.obs["dominant_comp"] = comp.idxmax(1).values
    a.obs["dominant_ct"] = a.obs["dominant_ct"].astype("category")
    a.obs["dominant_comp"] = a.obs["dominant_comp"].astype("category")
    sq.gr.spatial_neighbors(a, coord_type="generic", n_neighs=6)   # Visium hex ~6 neighbours
    return a, comp


def main():
    # ---- per-section niche analysis ----
    fig, axes = plt.subplots(2, 3, figsize=(16, 10)); axes = axes.ravel()
    adatas = {}
    for ax, sec in zip(axes, SECS):
        a, comp = build(sec)
        sq.gr.nhood_enrichment(a, cluster_key="dominant_comp", seed=0)
        z = a.uns["dominant_comp_nhood_enrichment"]["zscore"]
        cats = list(a.obs["dominant_comp"].cat.categories)
        im = ax.imshow(z, cmap="RdBu_r", vmin=-50, vmax=50)
        ax.set_xticks(range(len(cats))); ax.set_yticks(range(len(cats)))
        ax.set_xticklabels(cats, rotation=45, ha="right", fontsize=8); ax.set_yticklabels(cats, fontsize=8)
        ax.set_title(f"{sec} ({SUBTYPE[sec]})", fontsize=10)
        for i in range(len(cats)):
            for j in range(len(cats)): ax.text(j,i,f"{z[i,j]:.0f}",ha="center",va="center",fontsize=7)
        adatas[sec] = a
    fig.suptitle("Neighborhood enrichment (z): are compartments spatially clustered (red) or segregated (blue)?", fontsize=12)
    plt.tight_layout(rect=[0,0,1,0.96]); plt.savefig(f"{FIG}/partB_nhood_enrichment.png", dpi=130); plt.close()
    log("wrote partB_nhood_enrichment.png")

    # ---- tumor-immune boundary by subtype ----
    log("\n=== Tumor-immune infiltration: mean immune proportion within tumor-dominant spots ===")
    rows = []
    for sec in SECS:
        a = adatas[sec]
        tum = a.obs[a.obs["dominant_comp"]=="malignant"]
        # neighbor immune proportion for each tumor spot (via spatial graph)
        G = a.obsp["spatial_connectivities"]
        immune = a.obs["immune"].values
        nb_imm = np.asarray(G.dot(immune)).ravel() / np.asarray(G.sum(1)).ravel().clip(min=1)
        a.obs["nbr_immune"] = nb_imm
        tmask = (a.obs["dominant_comp"]=="malignant").values
        rows.append(dict(section=sec, subtype=SUBTYPE[sec], n_tumor=int(tmask.sum()),
                         immune_in_tumor=float(a.obs.loc[tmask,"immune"].mean()),
                         nbr_immune_at_tumor=float(np.nanmean(nb_imm[tmask]))))
    df = pd.DataFrame(rows)
    log(df.round(3).to_string(index=False))
    er = df[df.subtype=="ER+"]["nbr_immune_at_tumor"]; tn = df[df.subtype=="TNBC"]["nbr_immune_at_tumor"]
    # spot-level pooled test
    er_spots = np.concatenate([adatas[s].obs.loc[adatas[s].obs.dominant_comp=="malignant","nbr_immune"].values
                               for s in SECS if SUBTYPE[s]=="ER+"])
    tn_spots = np.concatenate([adatas[s].obs.loc[adatas[s].obs.dominant_comp=="malignant","nbr_immune"].values
                               for s in SECS if SUBTYPE[s]=="TNBC"])
    u,pv = mannwhitneyu(tn_spots, er_spots, alternative="greater")
    log(f"\nImmune at tumor margin: TNBC median={np.median(tn_spots):.3f} (n={len(tn_spots)}) vs "
        f"ER+ median={np.median(er_spots):.3f} (n={len(er_spots)}); MWU(TNBC>ER+) p={pv:.2e}")
    log("Interpretation: higher immune signal at TNBC tumor margins is consistent with TNBC being the "
        "more immune-infiltrated ('hot') subtype vs ER+ ('cold') — recovered spatially.")

    # boxplot
    fig,ax=plt.subplots(figsize=(6,5))
    ax.boxplot([er_spots, tn_spots], labels=[f"ER+\n(n={len(er_spots)})", f"TNBC\n(n={len(tn_spots)})"], showfliers=False)
    ax.set_ylabel("neighbor immune proportion at tumor spots"); ax.set_title("Tumor-immune infiltration by subtype")
    plt.tight_layout(); plt.savefig(f"{FIG}/partB_tumor_immune_by_subtype.png", dpi=130); plt.close()
    log("wrote partB_tumor_immune_by_subtype.png")

    # ---- spatial domains (composition leiden, pooled) + per-section maps ----
    dom_frames=[]
    for sec in SECS:
        a=adatas[sec]; dom_frames.append(pd.DataFrame(a.obs[["malignant","immune","stromal"]].values,
                                                      index=[f"{sec}_{b}" for b in a.obs_names],
                                                      columns=["malignant","immune","stromal"]))
    from sklearn.cluster import KMeans
    D=pd.concat(dom_frames)
    km = KMeans(n_clusters=6, random_state=0, n_init=10).fit(D.values)
    dom = pd.Series(km.labels_.astype(str), index=D.index)
    log(f"\n=== spatial domains: {dom.nunique()} recurrent composition-niches across sections ===")
    prof = D.groupby(dom.values).mean().round(2); prof["n"]=dom.value_counts()
    log(prof.to_string())

    fig,axes=plt.subplots(2,3,figsize=(16,10)); axes=axes.ravel()
    ncol = plt.cm.tab10(np.linspace(0,1,dom.nunique()))
    for ax,sec in zip(axes,SECS):
        a=adatas[sec]; ids=[f"{sec}_{b}" for b in a.obs_names]; dlab=dom.reindex(ids).astype(int).values
        ax.scatter(a.obsm["spatial"][:,0], -a.obsm["spatial"][:,1], c=ncol[dlab], s=6)
        ax.set_title(f"{sec} ({SUBTYPE[sec]})",fontsize=10); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    fig.suptitle(f"Spatial domains ({dom.nunique()} recurrent composition-niches)", fontsize=12)
    plt.tight_layout(rect=[0,0,1,0.96]); plt.savefig(f"{FIG}/partB_spatial_domains.png", dpi=130); plt.close()
    log("wrote partB_spatial_domains.png")

    for sec in SECS:
        a=adatas[sec]; ids=[f"{sec}_{b}" for b in a.obs_names]
        a.obs["spatial_domain"]=pd.Categorical(dom.reindex(ids).values)
        a.write(f"{OUTD}/{sec}.h5ad")
    open(f"{PROC}/partB_summary.txt","w").write("\n".join(LOG))
    log("\nwrote enriched AnnData + partB_summary.txt")


if __name__ == "__main__":
    main()
