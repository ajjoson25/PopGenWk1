import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

vcf_path = "AltaiNea.hg19_1000g.22.mod.vcf"
bin_size = 100_000

diff_per_bin = pd.Series(dtype="int64")
snp_matrix   = pd.DataFrame(dtype=np.int64)
chrom_label  = None

for chunk in pd.read_csv(
    vcf_path,
    sep="\t",
    comment="#",
    header=None,
    usecols=[0, 1, 3, 4],
    names=["#CHROM","POS","REF","ALT"],
    dtype={"#CHROM":"category", "POS":"int64", "REF":"object", "ALT":"object"},
    chunksize=1_000_000
):
    if chrom_label is None and not chunk.empty:
        chrom_label = str(chunk["#CHROM"].iloc[0])

    b = (chunk["POS"] // bin_size) * bin_size
    diff_counts = chunk.loc[chunk["ALT"] != "."].groupby(b).size()
    diff_per_bin = diff_per_bin.add(diff_counts, fill_value=0)

    snp = chunk[["REF","ALT"]].copy()
    snp.loc[snp["ALT"] == ".", "ALT"] = snp["REF"]
    mask = snp["REF"].str.len().eq(1) & snp["ALT"].str.len().eq(1)
    snp = snp.loc[mask]

    if not snp.empty:
        m = snp.groupby(["REF","ALT"]).size().unstack("ALT", fill_value=0)
        snp_matrix = snp_matrix.add(m, fill_value=0)

res = diff_per_bin.sort_index().astype("int64")
fig, ax = plt.subplots(figsize=(12, 4))
x = res.index.to_numpy(); w = bin_size
ax.bar(x, res.values, width=w, align="edge")
ax.set_xlabel("Genomic position (bp)")
ax.set_ylabel("# differences")
if chrom_label is not None:
    ax.set_title(f"Chr {chrom_label} — {bin_size//1000} kb bins")
ax.margins(x=0); plt.tight_layout()
fig.savefig("differences_per_100kb.png", dpi=300, bbox_inches="tight")
plt.close(fig)

with open("vcf_summary.txt", "w") as f:
    f.write(f"# Summary for {vcf_path}\n# Bin Size={bin_size}\n\n")
    f.write("## Differences per 100kb\n")
    res.rename("differences").to_csv(f, sep="\t", index_label="Position"); f.write("\n")
    f.write("## SNP 4x4 matrix (REF x ALT)\n")
    snp_matrix.astype("int64").to_csv(f, sep="\t", index_label="REF"); f.write("\n")

