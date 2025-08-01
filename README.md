# MicrobRNA: Microbial Profiling from RNA-Seq Data

Emerging evidence highlights the
potential role of the human microbiome in cancer development and progression, with
microbial transcripts in tumor tissues offering a previously unexplored source of tumor-associated antigens. 
Here, it is presented a computational pipeline designed to detect
microbial RNAs, specifically from Archaea, Bacteria, and Fungi, in RNA sequencing data from DLBCL samples. By adapting the MetaPhlAn profiling framework to
the characteristics of short-read paired RNA-seq data, we enable a more accurate identification
of microbial taxa. In addition, downstream filtering
steps were developed to enhance the robustness of microbial predictions, with configurable
thresholds that depend on user-defined input parameters. The main goal is to systematically screen for non-human reads, enabling
the identification of potential pathogen-derived sequences associated with cancer biology,
thus expanding the repertoire of strategies for treating tumors in a more precise and
personalized manner.

The pipeline step's consist in:
  1. Alignining RNA-seq reads to the reference genome and extract non-mapping reads
  2. Aligning unmapped reads to MetaPhlAn database
  3. BLAST filtering
  4. MetaPhlAn quantification step
  5. Coverage-based filtering and plotting


![pipeline overview](pipeline-overview.png)

---

## Requirements

- `Mamba`
- `STAR index` for the reference genome of choice


### To install Mamba
Go to [Miniforge installation](https://conda-forge.org/download/) and select the correct release based on your OS.
After the download is completed, run the installer in your MacOS/Linux terminal:
```bash
bash Miniforge3-Linux-x86_64.sh
```

Close and reopen your terminal to load conda, or run:
```bash
source ~/.bashrc   # or source ~/.zshrc depending on your shell
```

Finally, confirm it works:
```bash
conda --version
```

Once Miniconda is installed, install Mamba:
```bash
conda install mamba -n base -c conda-forge
```

### To index a reference genonme with STAR
After creating and activating the pipeline environment [(see instruction later)](#usage) choose the reference genome of choice. 
A complete database of genomes is [NCBI genomes website](https://www.ncbi.nlm.nih.gov/home/genomes/) where the user can both 
download directly from it or by using the [NCBI datasets command line tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)

For a better alignment of input reads to the index it is recommended to also specify the `gtf` file while creating the index

```bash
STAR --runMode genomeGenerate --genomeDir path/to/empty/directory --runThreadN INTEGER --genomeFastaFiles /path/to/reference.fasta --sjdbGTFfile /path/to/annotations.gtf
```

Since 

---

## Input

- A directory of paired-end FASTQ files, named with `_R1.fastq.gz` and `_R2.fastq.gz` suffixes
- A STAR genome index directory
- Path to the installed MetaPhlAn database

---

## Usage

It requires for a conda environment to be created. Use Mamba to create the environment
```bash
mamba env create -f environment.yml
mamba activate microbRNA
```

And the pipeline can be executed by running:
```bash
bash pipeline.sh --index-dir /path/to/genome/index --fastq-dir /path/to/fastq/files [options]
```

### Optional arugments
```text
-c, --coverage N	Coverage thresholds for marker positivity (default: 30 40)
-t, --threads N		Number of threads to use (default: 8)
-u, --unmapped-dir DIR	Output directory for unmapped reads (default: ./unmapped_reads)
-w, --working-dir DIR	Directory for intermediate/output files (default: ./)
-o, --outdir DIR	Final output directory (default: ./results)
-s, --seq-identity N	Sequence identity threshold(s) (default: 98 80)
-h, --help		Show usage help and exit
```

---

## Output
A binary infection matrix, where in each cell _(i, j)_ the value is `1` if the sample _i_ is infected by the species _j_ 

One binary matrix is output for each combination of the input parameters.

Also, a `plots` directory is created that contains as many subdirectories as the combinations of input parameteres, each containing the precomputed plots
on the corresponding output binary matrix. Each one of them is called `si{sequence-identity-threshold}-cov{coverage-threshold}`
