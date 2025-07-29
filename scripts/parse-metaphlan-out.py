import re
import os
import sys

# Sanity check
if len(sys.argv) > 1:
    working_dir = sys.argv[1]
    mpa_results_dir = sys.argv[2]
    si = sys.argv[3]
    final = int(sys.argv[4])
else:
    print("No argument provided.")
    sys.exit(0)


def parse_metaphlan_output(filename):
    """
    Legge il file di output di MetaPhlAn e restituisce:
      - total_reads: il numero totale di reads processate (estratto dalla riga "#... reads processed")
      - species_counts: un dizionario con chiave il nome della specie (la parte che inizia con "s__")
        e valore il count stimato (estratto dalla quinta colonna)
    """
    total_reads = None
    species_counts = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # Cerca la riga con il totale di reads processate
            if line.startswith('#'):
                m = re.search(r"#(\d+)\s+reads processed", line)
                if m:
                    total_reads = int(m.group(1))
                continue
            
            # Se la riga non è un commento, si assume che sia una riga della tabella
            # I campi sono separati da spazi o tab
            clade_name = line.split()[0]
            
            if len(clade_name.split('|')) != 7:
                continue
            
            # Considera solo le righe che contengono il livello specie (ovvero che contengono "|s__")
            if '|s__' in clade_name and 't__' not in clade_name:
                try:
                    # Il count stimato è il quinto campo (indice 4)
                    count = float(line.split()[4])
                except ValueError:
                    count = 0.0
                    
                # Estrae il nome della specie (si prende l'ultima parte che inizia con "s__")
                species_name = clade_name.split('|')[-1]
                species_name = species_name[3:]
              
                # if species_name.split('_')[1].startswith('SGB'):
                #     print(filename)
                
                species_counts[species_name] = count

    return total_reads, species_counts


# ====================================================================================================

### LOAD DATA FROM OUTPUT METAPHLAN TXT FILES 

# Dictionary that associates each sample to another dictionary, that associated each species found
# in that sample with the number of estimated counts per million
data = {}

samples = [f for f in os.listdir(mpa_results_dir) if os.path.isdir(os.path.join(mpa_results_dir, f))]

for sample in samples:
  
    sample_path = os.path.join(mpa_results_dir, sample)
    
    if os.path.isdir(sample_path):
        if final:
            file_path = os.path.join(sample_path, f'{sample}-final_{si}.txt')
        else:
            file_path = os.path.join(sample_path, f'profiled-{sample}.txt')
            
        
    if os.path.exists(file_path):
        # Retrieve total reads in the sample and read counts for each species detected
        total_reads, species_counts = parse_metaphlan_output(file_path)
        
        species_cpm = {}
        for species, count in species_counts.items():
            species_cpm[species] = (count / total_reads) * 1e6
  
        data[sample] = species_cpm


### GENERATE OUTPUT FILE

from collections import defaultdict

# Initialize dictionaries to count occurrences of each species
species_presence = defaultdict(set) 

# Iterate over data
for sample, species_dict in data.items():

    for species in species_dict.keys():
        species_presence[species].add(sample)


# Write to a file all the different species name so that in downstream analysis 
# (e.g. marker-coverage.sh) I can loop on all the species found in the dataset and have unique
# species names to a text file
if final:
    output_species_file = f'{working_dir}/postBLAST_species_file_{si}.tsv'
else:
    output_species_file = f'{working_dir}/preBLAST_species_file.tsv'


with open(output_species_file, 'w') as f:
    for species in sorted(species_presence.keys()):
        samples = sorted(species_presence[species])
        f.write(f"{species}\t{len(samples)}\t{','.join(samples)}\n")
        
        
'''        
df_plot = pd.DataFrame(data)

df_plot = df_plot.sort_values(by="Cell Lines Count", ascending=False)

# Step 3: Plot using Seaborn
plt.figure(figsize=(14, 6))
ax = sns.barplot(data=df_plot, x="Species", y="Cell Lines Count", hue="Organism", hue_order=all_categories, palette=palette_dict)

# Remove the box (spines)
for spine in ax.spines.values():
    spine.set_visible(False)
    
# Ensure only integer ticks on Y-axis
from matplotlib.ticker import MaxNLocator
ax.yaxis.set_major_locator(MaxNLocator(integer=True))

# Ensure tick labels are correctly aligned
ax.set_xticks(range(len(df_plot["Species"].unique())))  
ax.set_xticklabels(df_plot["Species"].unique(), rotation=60, ha="right")  

plt.ylabel("Number of Cell Lines")
plt.xlabel("Bacterial Species")
plt.title("Comparison of Bacterial Presence Across Cell Lines")

plt.tight_layout()
plt.savefig(f'{working-dir}/plots/bacterial_presence-preBlast.png')



# --- threshold: bacteria that are present in at least 3 cell lines
df_plot = df_plot.loc[df_plot["Cell Lines Count"] >= 4]

plt.figure(figsize=(14, 6))
ax = sns.barplot(data=df_plot, x="Species", y="Cell Lines Count", hue="Organism", hue_order=all_categories, palette=palette_dict)

# Remove the box (spines)
for spine in ax.spines.values():
    spine.set_visible(False)

# Ensure tick labels are correctly aligned
ax.set_xticks(range(len(df_plot["Species"].unique())))  
ax.set_xticklabels(df_plot["Species"].unique(), rotation=60, ha="right")  

plt.ylabel("Number of Cell Lines")
plt.xlabel("Bacterial Species")
plt.title("Comparison of Bacterial Presence Across Cell Lines")
plt.legend(title="Experiment", loc="upper right")

plt.tight_layout()

plt.savefig(f'{working-dir}/plots/bacterial_presence-preBlast_thr3.png')
'''


###################################################################################
### RIVERPLOT ###
'''
if final:
    import plotly.graph_objects as go
    
    # Step 1: Generate a list of unique node labels
    species_list = set()
    samples_list = set()
    
    for sample, species_dict in data.items():
        samples_list.add(sample)
        species_list.update(species_dict.keys())
    
    species_list = sorted(list(species_list))
    samples_list = sorted(list(samples_list))
    
    # Combine to create a full list of node names
    node_labels = samples_list + species_list
    
    # Create a mapping from label to index
    label_to_index = {label: idx for idx, label in enumerate(node_labels)}
    
    # Step 2: Build Sankey links (from sample to species)
    source = []
    target = []
    value = []
    
    for sample, species_dict in data.items():
        for species, cpm in species_dict.items():
            if cpm > 0:
                source.append(label_to_index[sample])
                target.append(label_to_index[species])
                value.append(cpm)  # You can also just use 1 if you want binary presence
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=node_labels
        ),
        link=dict(
            source=source,
            target=target,
            value=value
        ))])
        
    fig.update_layout(
        title_text="Microbial Species Across Cell Lines",
        font=dict(size=16),
        width=2000,
        height=1200
    )
    
    fig.update_layout(title_text="Microbial Species Across Cell Lines", font_size=10)
    #fig.write_html(f'{working_dir}/plots/riverplot.html')
    fig.write_html(f'{working_dir}/riverplot.html')
'''
