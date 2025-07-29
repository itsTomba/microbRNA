# MicrobRNA: Microbial Profiling from RNA-Seq Data

Emerging evidence highlights the
potential role of the human microbiome in cancer development and progression, with
microbial transcripts in tumor tissues offering a previously unexplored source of tumor-associated antigens. 
Here, it is presented a computational pipeline designed to detect
microbial nucleic acids, specifically from Archaea, Bacteria, and Fungi, in RNA sequencing data from DLBCL samples. By adapting the MetaPhlAn profiling framework to
the characteristics of RNA-seq data, this study enables a more accurate identification
of microbial taxa associated with lymphoma tissues. In addition, downstream filtering
steps were developed to enhance the robustness of microbial predictions, with configurable
thresholds that depend on user-defined input parameters. The main goal is to systematically screen for non-human reads, enabling
the identification of potential pathogen-derived sequences associated with cancer biology,
thus expanding the repertoire of strategies for treating tumors in a more precise and
personalized manner.

It is designed to:
1. Align RNA-seq reads to the human genome and extract non-mapping reads
2. Align reads to MetaPhlAn database
3. BLAST filtering
4. MetaPhlAn quantification step
5. Coverage-based filtering and plotting

![pipeline overview](pipeline-overview.PNG)

## Requirements

- Reference genome indexed
- `STAR`
- `samtools`
- `Conda`


## Input

- A directory of FASTQ paired-end files, named with `_R1.fastq.gz` and `_R2.fastq.gz` suffixes
- A STAR genome index directory
- A MetaPhlAn database installation


## Usage

It requires for a conda environment to be created:
```bash
conda env create -f environment.yml
conda activate microbRNA
```

And the pipelince can be executed by running:
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

## Output
A binary infection matrix, where in each cell _(i, j)_ the value is `1` if the species _i_ is infected by the species _j_ 

One binary matrix is output for each combination of the input parameters.

Also, a `plots` directory is created that contains as many subdirectories as the combinations of input parameteres, each containing the precomputed plots
on the corresponding output binary matrix. Each one of them is called `si{sequence-identity-threshold}-cov{coverage-threshold}`.
