################################################################################
### NEW PLOT TO VISUALIZE A LOT OF SAMPLES
################################################################################

import os
import sys
import pandas as pd
from collections import defaultdict


# Usage check
if len(sys.argv) > 1:
    working_dir = sys.argv[1]
    si = sys.argv[2]
    coverage = sys.argv[3]
else:
    print(f"Usage: {sys.argv[0]} <working_dir> <s.i. BLAST> <coverage>")
    sys.exit(1)

plots_dir = os.path.join(working_dir, "plots", f"si{si}-cov{coverage}")
os.makedirs(plots_dir, exist_ok=True)


from ete3 import NCBITaxa

ncbi = NCBITaxa()

def organism_type(species_name):
    species_name = species_name.replace("[", "").replace("]", "")
    
    if 'GGB' in species_name:
        return "Unknown"
      
    parts = species_name.split("_")
    
    if 'SGB' in species_name or len(parts) < 2 or parts[1] == 'sp':
        search_name = parts[0]
    else:
        search_name = " ".join(parts[:2])
    
    try:
        taxid = ncbi.get_name_translator([search_name])[search_name][0]
    except KeyError:
        print("Species not found with 2-part name:", search_name)
        try:
            other_name = search_name.replace("us", "oides")
            taxid = ncbi.get_name_translator([other_name])[other_name][0]
        except KeyError:
            print("Also not found with other name:", other_name)
            return "Unknown"
      
    lineage = ncbi.get_lineage(taxid)
    
    if 2 in lineage:
        return "Bacteria"
    elif 4751 in lineage:
        return "Fungi"
    elif 2157 in lineage:
        return "Archea"
    else:
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        for tid in lineage:
            if ranks[tid] == "kingdom":
                return names[tid]
        return names.get(lineage[1], "Unknown")


############### Parse the filtered_infections.tsv file

'''
# Read the file
file_path = f"{working_dir}/presence_absence_matrix_{si}_{int(coverage)}.csv"
df = pd.read_csv(file_path, index_col=0)
'''

'''
df_long = pd.read_csv(f"{working_dir}/filtered_infections_file_{si}_{coverage}.tsv", sep="\t")  # Adjust filename/path & sep if needed

# Group by species: get number of unique samples and list of samples
summary = df_long.groupby('species').agg(
    num_samples=('sample', 'nunique'),
    samples_list=('sample', lambda x: ",".join(sorted(x.unique())))
).reset_index()
summary.to_csv(f"{working_dir}/converted_filtered_infections.tsv", sep="\t", index=False, header=False)
'''


with open(f"{working_dir}/filtered_infections_file_{si}_{coverage}.tsv") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        species, count, sample_list = parts
        samples = sample_list.split(",")
        data[species].update(samples)
        all_samples.update(samples)

# Create binary presence/absence matrix
all_samples = sorted(all_samples)
df = pd.DataFrame(0, index=data.keys(), columns=all_samples)

for species, samples in data.items():
    for sample in samples:
        df.loc[species, sample] = 1


############### Histogram of per-sample species richness

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator   # to ensure integer values as x-label ticks
import seaborn as sns

# Calculate species richness per sample
species_counts = df.sum(axis=0)

# Identify samples with the highest richness
max_species = int(species_counts.max())

# Plot with seaborn
plt.figure(figsize=(8, 6))
sns.histplot(species_counts, bins=range(1, max_species + 2), color="royalblue", edgecolor="black")
sns.despine()

plt.xlabel("Number of Microbial Species")
plt.ylabel("Number of Samples")
plt.title("Species number per Sample")

# Set integer ticks on the x-axis
ax = plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.tight_layout()
plt.savefig(f"{plots_dir}/per-sample-species-richness_{si}_{coverage}.png", dpi=300)


############### Barplot of most widespread species

# Prepare the data
sample_counts = df.sum(axis=1).sort_values(ascending=False)
top_species = sample_counts.head(20)

top_species_df = top_species.reset_index()
top_species_df.columns = ["Species", "SampleCount"]

# Assign kingdom to each species
top_species_df["Kingdom"] = top_species_df["Species"].apply(organism_type)

# Plot
plt.figure(figsize=(10, 6))
sns.barplot(data=top_species_df, x="SampleCount", y="Species", hue="Kingdom", palette="Set2")
sns.despine()
plt.xlabel("Number of Samples")
#plt.ylabel("Species")
plt.title("Top 20 Most Prevalent Microbial Species")
plt.legend(title="Kingdom")
plt.tight_layout()
plt.savefig(f"{plots_dir}/most-common-species_{si}_{coverage}.png")


############### Presence/absence heatmap with hierarchical clustering

import numpy as np

# I have a species that is present in all samples, so its row on the binary 
# matrix will have all ones. This makes the clustering algorithm throw an error
# because that row is constant. So for the scope of this visualization, I have to 
# remove that row (that species)
  # print((df.nunique(axis=1) == 1).sum(), "rows have constant values")
  # print((df.nunique(axis=0) == 1).sum(), "columns have constant values")
# Filter out rows that have constant values across samples
df_filtered = df.loc[df.nunique(axis=1) > 1]

g = sns.clustermap(df_filtered, cmap="Greens", metric="euclidean", method="average", standard_scale=0, figsize=(15, 12), cbar_kws={"label": "Presence"}, row_cluster=False )
plt.savefig(f"{plots_dir}/heatmap_{si}_{coverage}.png")


############### Species discovery rate

 # Clears the current figure so that we can do another plot
plt.clf()
plt.figure(figsize=(10, 6))

discovery_curve = []
seen_species = set()
for sample in df.columns:
    new_species = df[sample][df[sample] == 1].index
    seen_species.update(new_species)
    discovery_curve.append(len(seen_species))

# Prepare data for Seaborn
curve_df = pd.DataFrame({
    "Number of Samples": range(1, len(df.columns) + 1),
    "Cumulative Unique Species": discovery_curve
})

# Plot using Seaborn
sns.lineplot(data=curve_df, x="Number of Samples", y="Cumulative Unique Species", color="royalblue")
sns.despine()
plt.xlabel("Number of Samples")
plt.ylabel("Cumulative Unique Species")
plt.title("Species Discovery Curve")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{plots_dir}/discovery-curve_{si}_{coverage}.png")


############### Pie chart

# Step 1: Count number of samples per species
df_prevalence = pd.DataFrame({
    "species": df.index,
    "num_samples": df.sum(axis=1)
})

total_samples = 38

# Thresholds
rare_threshold = int(total_samples * 0.05)          # ≤5% → rare
intermediate_min = rare_threshold + 1
intermediate_max = int(total_samples * 0.5)         # >5% and ≤50% → intermediate
high_threshold = intermediate_max                   # >50% → plotted individually

# Step 1: Count number of samples per species
df_prevalence = pd.DataFrame({
    "species": df.index,
    "num_samples": df.sum(axis=1)
})

# Step 2: Categorize
rare_df = df_prevalence[df_prevalence["num_samples"] <= rare_threshold]
intermediate_df = df_prevalence[
    (df_prevalence["num_samples"] >= intermediate_min) &
    (df_prevalence["num_samples"] <= intermediate_max)
]
# Each high-prevalence species plotted separately
individual_df = df_prevalence[df_prevalence["num_samples"] > high_threshold]

# Step 3: Build pie chart data
labels = individual_df["species"].tolist()
sizes = individual_df["num_samples"].tolist()

if not intermediate_df.empty:
    labels.append(f"Intermediate species ({intermediate_min}–{intermediate_max} samples)")
    sizes.append(intermediate_df["num_samples"].sum())

if not rare_df.empty:
    labels.append(f"Rare species (≤{rare_threshold} samples)")
    sizes.append(rare_df["num_samples"].sum())

# Step 5: Plot

# Base palette for individual species
base_colors = plt.cm.tab20.colors

# Fixed colors for categories
rare_color = "lightcoral"
intermediate_color = "skyblue"

# Assign colors to each slice
colors = []
idx = 0
for lbl in labels:
    if "Rare species" in lbl:
        colors.append(rare_color)
    elif "Intermediate species" in lbl:
        colors.append(intermediate_color)
    else:
        colors.append(base_colors[idx % len(base_colors)])
        idx += 1

plt.figure(figsize=(20, 12))


# Optional explode last two slices
explode = [0.1 if "Intermediate" in label or "Rare" in label else 0 for label in labels]
    
wedges, texts, autotexts = plt.pie(
    sizes,
    labels=labels,
    autopct='%1.1f%%',
    startangle=140,
    wedgeprops=dict(width=0.4),
    colors=colors[:len(labels)],
    explode=explode,
    textprops={'fontsize': 16}
)

# Adjust text size
for t in texts:
    t.set_fontsize(18)
for at in autotexts:
    at.set_fontsize(16)

plt.tight_layout()
plt.savefig(f"{plots_dir}/species_prevalence_pie_{si}_{coverage}.png", dpi=300)


### Version with no label so i can write them bigger
plt.figure(figsize=(20, 12))
wedges, texts, autotexts = plt.pie(
    sizes,
    #labels=labels,
    autopct='%1.1f%%',
    startangle=140,
    wedgeprops=dict(width=0.4),
    colors=colors[:len(labels)],
    explode=explode,
    textprops={'fontsize': 16}
)

# Adjust text size
for t in texts:
    t.set_fontsize(0)
for at in autotexts:
    at.set_fontsize(0)

plt.tight_layout()
plt.savefig(f"{plots_dir}/species_prevalence_pie_{si}_{coverage}_noLabel.png", dpi=300)


'''
############### Co-occurrence circle plot

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from itertools import combinations
from collections import defaultdict
import math

all_samples = sorted(all_samples)
df = pd.DataFrame(0, index=data.keys(), columns=all_samples)
for species, samples in data.items():
    df.loc[species, list(samples)] = 1

# Step 3: Build the co-occurrence network
G = nx.Graph()

# Only consider species that appear in the dataset
G.add_nodes_from(df.index)

# Add edges based on co-occurrence across samples
for sample in df.columns:
    species_present = df[df[sample] == 1].index.tolist()
    for sp1, sp2 in combinations(species_present, 2):
        if G.has_edge(sp1, sp2):
            G[sp1][sp2]['weight'] += 1
        else:
            G.add_edge(sp1, sp2, weight=1)

# Step 4: Filter edges by co-occurrence threshold (≥ 20 samples)
threshold = 100
strong_edges = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] >= threshold]
filtered_nodes = set(u for edge in strong_edges for u in edge)

# Create subgraph
G_filtered = G.edge_subgraph(strong_edges).copy()
G_filtered.add_nodes_from(filtered_nodes)  # ensure isolated nodes remain

# Step 5: Draw the circular plot
plt.figure(figsize=(14, 14))
pos = nx.circular_layout(G_filtered)

# Draw nodes
nx.draw_networkx_nodes(G_filtered, pos, node_size=700, node_color='lightgreen', edgecolors='black')

# Draw edges with width proportional to weight
edges = G_filtered.edges(data=True)
nx.draw_networkx_edges(G_filtered, pos, edgelist=edges, width=[d['weight'] / 5 for (_, _, d) in edges], alpha=0.6)

# Draw labels
# Draw outer labels (no rotation)
for node, (x, y) in pos.items():
    angle = math.atan2(y, x)
    label_distance = 1.1
    label_x = label_distance * math.cos(angle)
    label_y = label_distance * math.sin(angle)
    plt.text(
        label_x,
        label_y,
        node,
        fontsize=16,
        ha='center',
        va='center'
    )

plt.title(f"Microbial Co-occurrence Network (≥ {threshold} Shared Samples)", fontsize=16)
plt.axis('off')
plt.tight_layout()
plt.savefig(f"{working_dir}/plots/species_cooccurrence_network_{si}.png", dpi=300)
'''

############### Binary co-occurence matrix

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
import matplotlib.ticker as ticker

# Compute co-occurrence matrix
cooccurrence_matrix = df.dot(df.T)

# Remove diagonal (self-cooccurrence)
np.fill_diagonal(cooccurrence_matrix.values, 0)

# Select top 50 species by total co-occurrence count
top_species = cooccurrence_matrix.sum(axis=1).sort_values(ascending=False).head(15).index
reduced_matrix = cooccurrence_matrix.loc[top_species, top_species]

# Create mask for diagonal
mask = np.eye(len(reduced_matrix), dtype=bool)

# Plot the heatmap
plt.figure(figsize=(14, 12))
ax = sns.heatmap(
    reduced_matrix,
    cmap="Greens",
    linewidths=0.5,
    linecolor="gray",
    square=True,
    cbar_kws={"label": "Co-occurrence Count"},
    mask=mask,
    annot=False  # You can set to True to display counts in cells
)

# Black out the diagonal with a second plot on top (optional visual clarity)
#for i in range(len(species)):
#    plt.gca().add_patch(plt.Rectangle((i, i), 1, 1, color='black'))


# Set font size for labels
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=14)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=14)

# Force colorbar ticks to be integers
colorbar = ax.collections[0].colorbar
colorbar.ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

plt.tight_layout()
plt.savefig(f"{plots_dir}/top15_cooccurrence_matrix_{si}_{coverage}.png", dpi=300)


############### Venn diagram by kingdom

import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Step 1: Load presence/absence matrix
matrix_df = pd.read_csv(f"{working_dir}/presence_absence_matrix_{si}_{int(coverage)}.csv", index_col=0)

# Step 2: Convert to long format (species, sample, presence)
long_df = matrix_df.reset_index().melt(id_vars="species", var_name="sample", value_name="present")
long_df = long_df[long_df["present"] == 1]  # Keep only present cases

# Step 3: Assign kingdoms
long_df["kingdom"] = long_df["species"].apply(organism_type)

# Step 4: Build sets of samples by kingdom
bacteria_samples = set(long_df[long_df["kingdom"] == "Bacteria"]["sample"])
fungi_samples = set(long_df[long_df["kingdom"] == "Fungi"]["sample"])
unknown_samples = set(long_df[long_df["kingdom"] == "Unknown"]["sample"])

# Step 5: Plot Venn diagram
plt.figure(figsize=(8, 6))
venn3(
    [bacteria_samples, fungi_samples, unknown_samples],
    set_labels=("Bacteria", "Fungi", "Unknown")
)
plt.title("Overlap of Infected Samples by Kingdom")
plt.tight_layout()
plt.savefig(f"{plots_dir}/venn-by-kingdoms_{si}_{coverage}.png")

