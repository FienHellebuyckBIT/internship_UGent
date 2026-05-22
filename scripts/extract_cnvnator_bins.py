import uproot
import pandas as pd
import numpy as np
import sys

#Import arguments from bash script
ROOT_FILE = sys.argv[1]
BIN_SIZE = int(sys.argv[2])
tsv_file = sys.argv[3]

f = uproot.open(ROOT_FILE)

all_data = []

chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

for chrom in chromosomes:

    hist_name = f"bin_{BIN_SIZE}/his_rd_p_{chrom}_{BIN_SIZE}"

    print(f"Reading: {hist_name}")

    try:

        hist = f[hist_name]

        values = hist.values()
        edges = hist.axis().edges()

        starts = edges[:-1].astype(int)
        ends = edges[1:].astype(int)

        df = pd.DataFrame({
            "chr": chrom,
            "start": starts,
            "end": ends,
            "rd": values
        })

        all_data.append(df)

    except Exception as e:

        print(f"SKIP {chrom}: {e}")

final_df = pd.concat(all_data, ignore_index=True)

# Remove invalid bins
final_df = final_df[
    np.isfinite(final_df["rd"])
]

final_df = final_df[
    final_df["rd"] > 0
]

# Convert to log2 ratio
genome_mean = final_df["rd"].mean()
final_df["log2ratio"] = np.log2(final_df["rd"] / genome_mean)


# Midpoint
final_df["position"] = (
    final_df["start"] +
    final_df["end"]
) / 2

final_df.to_csv(
    tsv_file,
    sep="\t",
    index=False
)

print()
print(final_df.head())
print()
print(f"TOTAL BINS: {len(final_df):,}")
print("DONE")
