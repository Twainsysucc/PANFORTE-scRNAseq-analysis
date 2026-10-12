# 01_QC.py — OmicVerse / Scanpy QC

from pathlib import Path
import anndata as ad
import omicverse as ov
import scanpy as sc

data_dir = Path("data")
sample_dirs = sorted(path for path in data_dir.iterdir() if path.is_dir())

adatas = []
for sample_dir in sample_dirs:
    sample = sc.read_10x_mtx(
        sample_dir,
        var_names="gene_symbols",
        cache=True
    )
    sample.obs["sample_id"] = sample_dir.name
    adatas.append(sample)

adata = ad.concat(adatas, merge="same")
adata.var_names_make_unique()
adata.obs_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    percent_top=[20],
    log1p=True,
    inplace=True
)
print(f"Cells before QC: {adata.n_obs}")

adata_filter = ov.pp.qc(
    adata,
    mode="mads",
    nmads=5,
    min_cells=3,
    min_genes=200,
    tresh={
        "mito_perc": 0.7,
        "nUMIs": 400,
        "detected_genes": 400
    }
)
print(f"Cells after QC: {adata_filter.n_obs}")



