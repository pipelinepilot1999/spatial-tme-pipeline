"""
06_figures.py — Hero figure: each Visium section colored by deconvolved composition.
RGB encoding per spot: R=malignant, G=immune, B=stromal (proportions sum to 1, so the
color IS the composition). Paired with the pathologist map for visual validation.
Outputs: figures/hero_composition.png, figures/composition_vs_pathology.png
"""
import os, numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PROC = "/home/ubuntu/spatial-tme-pipeline/data/processed"
FIG  = "/home/ubuntu/spatial-tme-pipeline/figures"; os.makedirs(FIG, exist_ok=True)
SECS = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]
SUBTYPE = {"1142243F":"TNBC","1160920F":"TNBC","CID4290":"ER+","CID4465":"TNBC","CID44971":"TNBC","CID4535":"ER+"}
COMPART = {"ER+ tumor":"malignant","TNBC tumor":"malignant","T cell":"immune",
           "Myeloid":"immune","Plasma cell":"immune","Fibroblast":"stromal","Endothelial":"stromal"}


def comp_props(sec):
    p = pd.read_csv(f"{PROC}/rctd_output/{sec}_rctd_proportions.csv", index_col=0)
    p.index = p.index.astype(str)
    out = pd.DataFrame(index=p.index)
    for comp in ["malignant","immune","stromal"]:
        out[comp] = p[[c for c in p.columns if COMPART.get(c)==comp]].sum(1)
    return out.div(out.sum(1).replace(0,np.nan), axis=0).fillna(0)


def coords(sec):
    m = pd.read_csv(f"{PROC}/rctd_input/{sec}/spot_meta.csv", index_col=0)
    m.index = m.index.astype(str); return m


def hero():
    fig, axes = plt.subplots(2, 3, figsize=(15, 10)); axes = axes.ravel()
    for ax, sec in zip(axes, SECS):
        c = comp_props(sec); m = coords(sec).loc[c.index]
        rgb = c[["malignant","immune","stromal"]].values.clip(0,1)
        ax.scatter(m["x"], -m["y"], c=rgb, s=6)
        ax.set_title(f"{sec} ({SUBTYPE[sec]})", fontsize=11); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    fig.suptitle("Breast cancer Visium — deconvolved composition (RCTD)\n"
                 "R = malignant   G = immune   B = stromal", fontsize=14)
    fig.text(0.5, 0.02, "Each spot colored by its cell-type composition; the spatial 'where' "
             "bulk & dissociated scRNA-seq cannot resolve.", ha="center", fontsize=9, style="italic")
    plt.tight_layout(rect=[0,0.03,1,0.95]); plt.savefig(f"{FIG}/hero_composition.png", dpi=140); plt.close()
    print("wrote hero_composition.png")


def vs_pathology():
    fig, axes = plt.subplots(2, 6, figsize=(22, 8))
    for j, sec in enumerate(SECS):
        c = comp_props(sec); m = coords(sec).loc[c.index]
        rgb = c[["malignant","immune","stromal"]].values.clip(0,1)
        axes[0,j].scatter(m["x"], -m["y"], c=rgb, s=5); axes[0,j].set_title(f"{sec}\n{SUBTYPE[sec]}", fontsize=9)
        # pathology: malignant-containing = red, stroma = green, immune = blue, other = grey
        lab = m["pathology"].astype(str).str.lower()
        col = np.full((len(m),3), 0.7)
        col[lab.str.contains("stroma|adipose").values] = [0,0.7,0]
        col[lab.str.contains("lymphocyte|tls").values] = [0,0,0.9]
        col[lab.str.contains("invasive cancer|dcis|cancer trapped").values] = [0.9,0,0]
        axes[1,j].scatter(m["x"], -m["y"], c=col, s=5); axes[1,j].set_title("pathology", fontsize=9)
        for r in (0,1): axes[r,j].set_xticks([]); axes[r,j].set_yticks([]); axes[r,j].set_aspect("equal")
    axes[0,0].set_ylabel("RCTD composition", fontsize=11); axes[1,0].set_ylabel("pathologist", fontsize=11)
    plt.suptitle("Deconvolved composition vs pathologist annotation (red=cancer, green=stroma, blue=immune)", fontsize=13)
    plt.tight_layout(); plt.savefig(f"{FIG}/composition_vs_pathology.png", dpi=130); plt.close()
    print("wrote composition_vs_pathology.png")


if __name__ == "__main__":
    hero(); vs_pathology()
