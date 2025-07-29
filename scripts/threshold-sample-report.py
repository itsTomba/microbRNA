################################################################################
##  Looks into the final predictions from final_species_file.tsv, filters them 
##  w.r.t. the coverage of their markers with a threshold given by the user, 
##  and outputs a summary of the infections (most widespread species, number 
##  of infected samples, ...)
################################################################################

import os
import sys
import pandas as pd

# Usage check
if len(sys.argv) > 1:
    working_dir = sys.argv[1]
    mpa_results_dir = sys.argv[2]
    si = sys.argv[3]
    coverage = float(sys.argv[4])
else:
    print(f"Usage: {sys.argv[0]} <working_dir> <mpa_results_dir> <s.i. BLAST> <coverage>")
    sys.exit(1)

import pandas as pd

species_file = f"{working_dir}/postBLAST_species_file_{si}.tsv"
species_df = pd.read_csv(species_file, sep="\t", header=None, names=["species", "count", "samples"])

# Build sample → species mapping
sample_to_species = {}
for _, row in species_df.iterrows():
    for sample in row["samples"].split(","):
        sample_to_species.setdefault(sample.strip(), []).append(row["species"])

def simplify_rname(rname):
    parts = rname.split('|')
    if len(parts) == 3 and parts[1].startswith('GC'):
        return f"{parts[0][0]}.{parts[0].split('_')[1]} {parts[2]}"
    elif len(parts) == 1 and '__' in rname:
        return "-".join(rname.split("__")[1].split("-")[:2])
    return rname

# For collecting confirmed infections
high_confidence_infections = []

for sample in sorted(sample_to_species.keys()):
    sample_path = os.path.join(mpa_results_dir, sample)
    coverage_file = f"{sample_path}/new-coverage/{sample}-coverage-postBLAST_{si}.tsv"

    if not os.path.isfile(coverage_file):
        print(f"[!] Missing coverage file for sample {sample}, skipping.")
        continue

    print(f"[*] Processing {sample}...")
    try:
        cov_df = pd.read_csv(coverage_file, sep="\t", usecols=["#rname", "coverage"])
    except Exception as e:
        print(f"[!] Failed to read coverage file for {sample}: {e}")
        continue

    cov_df["#rname"] = cov_df["#rname"].apply(simplify_rname)

    confirmed_species = []

    for species in sample_to_species.get(sample, []):
        marker_file = f"{working_dir}/species-markers/markers-{species}.txt"
        if not os.path.isfile(marker_file):
            continue

        with open(marker_file) as f:
            markers = set(line.strip() for line in f if line.strip())
        simplified_markers = set(simplify_rname(m) for m in markers)

        matched = cov_df[cov_df["#rname"].isin(simplified_markers)]
        if matched.empty:
            continue

        if (matched["coverage"] > coverage).any():
            confirmed_species.append(species)
            high_confidence_infections.append((species, sample))

    print(f"  ↳ {len(confirmed_species)} confident species")

# Build species → list of samples mapping
species_to_samples = {}
for species, sample in high_confidence_infections:
    species_to_samples.setdefault(species, set()).add(sample)


# Create binary presence/absence matrix
exploded_rows = []
for species, samples in species_to_samples.items():
    for sample in samples:
        exploded_rows.append({"species": species, "sample": sample})

presence_df = pd.DataFrame(exploded_rows)
presence_df["present"] = 1

matrix_df = (
    presence_df.pivot_table(index="species", columns="sample", values="present", fill_value=0)
    .sort_index(axis=1)
)

output_csv = os.path.join(working_dir, f"presence_absence_matrix_{si}_{int(coverage)}.csv")
matrix_df.to_csv(output_csv)

print(f"\n[✓] High-confidence presence/absence matrix saved to {output_csv}")



'''
# Save in the original format
output_file = os.path.join(working_dir, f"filtered_infections_file_{si}_{int(coverage)}.tsv")
with open(output_file, "w") as out:
    for species, samples in species_to_samples.items():
        sample_list = ",".join(sorted(samples))
        out.write(f"{species}\t{len(samples)}\t{sample_list}\n")

print(f"\n[✓] High-confidence infection report saved to {output_file}")

# Summary statistics
confirmed_df = pd.DataFrame(high_confidence_infections, columns=["species", "sample"])
n_infected_samples = confirmed_df["sample"].nunique()
n_species_found = confirmed_df["species"].nunique()

print("\n========== SUMMARY ==========")
print(f"Threshold for detection: ≥{coverage}% coverage on any marker")
print(f"Total high-confidence infections: {len(confirmed_df)}")
print(f"Number of infected samples: {n_infected_samples}")
print(f"Number of species detected: {n_species_found}")
print("=============================\n")
'''

# Optional: save a species × sample matrix
# matrix = confirmed_df.pivot_table(index="sample", columns="species", aggfunc=lambda x: 1, fill_value=0)
# matrix_file = os.path.join(working_dir, "species_sample_matrix.tsv")
# matrix.to_csv(matrix_file, sep="\t")
# print(f"[✓] Species × Sample matrix saved to {matrix_file}")
