#!/bin/bash

echo
echo "*************************************************************************"
echo "***************** MICROBIAL PROFILING FROM RNA-SEQ DATA *****************"
echo "*************************************************************************"
echo

################################################################################
# STEP 0 : CHECK PARAMETERS

print_help() {
cat <<EOF
Usage: $0 [OPTIONS]

Required arguments:
  -n, --star-index DIR      Path to the STAR reference genome index directory.
  -i, --fastq-dir DIR       Directory containing input FASTQ files.
  -d, --mpa-dir DIR         Path to the MetaPhlAn database.
  -b, --blast-dir DIR       Path to the BLAST database.

Optional arguments:
  -c, --coverage N          Percentage of coverage to consider a marker as positive. [default: (30 40)]
  -l, --blast-index STR     Name of BLAST database. [default: core_nt]
  -x, --mpa-index STR       Version of the MetaPhlAn database. [default: mpa_vJan25_CHOCOPhlAnSGB_202503 (latest)]
  -t, --threads N           Number of threads to use. [default: 8]
  -u, --unmapped-dir DIR    Directory where to store unmapped reads on reference genome. [default: ./unmapped_reads]
  -w, --working-dir DIR     Directory where intermediate files and output files will be store. [default: .]
  -o, --outdir DIR          Directory storing output files. [defaults: ./results]
  -s, --seq-identity N      Sequence identity threshold for considering a read human if blasted on a sequence. [default: (98 80) ]
  -h, --help                Show this help message and exit.

Example:
  $0 --star-index /path/to/index --fastq-dir /path/to/fastq/files --mpa-dir /path/to/metaphlan/database --blast-dir /path/to/blast/db/ --threads 8

EOF

exit 0
}

# Default values
star_index=""
fastq_dir=""
mpa_db=""
mpa_index="mpa_vJan25_CHOCOPhlAnSGB_202503"
blast_dir=""
blast_index="core_nt"
threads=8
unmapped_dir="$PWD/unmapped_reads"
wd="$PWD/"
plots_dir="$PWD/plots"
si_thr=("98" "80")
outdir="$PWD/results"
coverage=("30" "40")


# Manual parsing loop
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -i|--fastq-dir) fastq_dir="$2"; shift ;;
    -d|--mpa-dir) mpa_db="$2"; shift ;;
    -x|--mpa-index) mpa_index="$2"; shift ;;
    -b|--blast-dir) blast_dir="$2"; shift ;;
    -l|--blast-index) blast_index="$2"; shift ;;
    -n|--star-index) star_index="$2"; shift ;;
    -t|--threads) threads="$2"; shift ;;
    -u|--unmapped-dir) unmapped_dir="$2"; shift ;;
    -w|--working-dir) wd="$2"; shift ;; 
    -s|--seq-identity) 
            shift
            si_thr=()
            while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
              si_thr+=("$1")
              shift
            done
            continue ;;
    -c|--coverage)
            shift
            coverage=()
            while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
              coverage+=("$1")
              shift
            done
            continue ;;
    -o|--outdir) outdir="$2"; shift ;;
    -h|--help) print_help ; shift ;;
    *) echo -e "Unknown parameter passed: $1\n"; exit 1 ;;
  esac
  shift
done


# Check required directories
if [[ -z "$star_index" || -z "$fastq_dir" ]]; then
  echo -e "[FATAL] --index-dir and --fastq-dir are required\n"
  exit 1
fi


# Print the choose parameters:
echo "Chosen parameters:"
echo -e "  Genome index directory: $star_index"
echo -e "  Input directory with fastq files: $fastq_dir"
echo -e "  MetaPhlAn database directory: $mpa_db"
echo -e "  MetaPhlAn database version: $mpa_index"
echo -e "  BLAST database directory: $blast_dir"
echo -e "  BLAST database version: $blast_index"
echo -e "  Threads: $threads"
echo -e "  Directory storing unmapped_reads: $unmapped_dir"
echo -e "  Working directory: $wd"
echo -e "  Output directory: $outdir"
echo -e "  Sequence identity threshold(s):\n$(for si in "${si_thr[@]}"; do echo -e "\t- $si"; done)"
echo -e "  Coverage threshold to consider a marker 'positive':\n$(for cov in "${coverage[@]}"; do echo -e "\t- $cov"; done)"


# Check the existance of the input directory:
if [ ! -d "$fastq_dir" ]; then
  echo -e "[FATAL] Directory '$fastq_dir' does not exist.\n"
  exit 1
fi

# Check the existance of the genome index directory:
if [ ! -d "$star_index" ]; then
  echo -e "[FATAL] Directory '$star_index' does not exist."
  echo "[?] Generate one with the command:"
  echo "[?]   STAR --runMode genomeGenerate --genomeDir path/to/new/index --runThreadN 16 --genomeFastaFiles /path/to/reference/fasta --sjdbGTFfile /path/to/annotations.gtf"
  exit 1
fi

# Check the existance of unmapped reads directory, if not existing then asks the user if he wants to create it:
if [ ! -d "$unmapped_dir" ]; then
  echo -e "[!] Directory '$unmapped_dir' does not exist"
  read -p "[?] Do you want to create it now? (y -> continue /n -> exits ): " answer
  case "$answer" in
    [Yy]* )
      echo -e "[+] Creating '$unmapped_dir'...\n"
      mkdir -p "$unmapped_dir" ;;
    * )
      echo -e "[x] Exiting script.\n"
      exit 1 ;;
  esac
fi

# Check the existance of working directory, if not existing then asks the user if he wants to create it:
if [ ! -d "$wd" ]; then
  echo -e "[!] Directory '$wd' does not exist"
  read -p "[?] Do you want to create it now? (y/n): " answer
  case "$answer" in
    [Yy]* )
      echo -e "[+] Creating '$wd'...\n"
      mkdir -p "$wd" ;;
    * )
      echo -e "[x] Exiting script.\n"
      exit 1 ;;
  esac
fi


#############################################################################################################################################################################
# STEP 1 : STAR ALIGNMENT TO RETRIEVE NON-MAPPING READS

echo
echo "####################################"
echo "STEP 1: RETRIEVAL OF UNMAPPED READS"

for file in $(ls "$fastq_dir"/*_R1.fastq.gz );
do

  # Start counting time
  start=$(date +%s)
    
  sample="$( basename "$file" _R1.fastq.gz )"
  
    
  # Check if it is necessary to skip the alignment on current sample
  if [[ -f "${unmapped_dir}/$sample/${sample}_R1.Unmapped.out.gz" && -f "${unmapped_dir}/$sample/${sample}_R2.Unmapped.out.gz" ]]; then
      echo -e "[✓] Skipping $sample"
      continue
  fi

  
  echo
  echo "----- ${sample} -----"
  
  read1="$file"
  read2="$fastq_dir/${sample}_R2.fastq.gz"
  # echo "read1: $read1"
  # echo "read2: $read2"
  
  mkdir -p "${unmapped_dir}/${sample}"
  
  if [[ ! -f "${unmapped_dir}/$sample/${sample}.Aligned.sortedByCoord.out.bam" ]]; then
      echo "[ ] STAR alignment"
      STAR --runMode alignReads \
        --readFilesIn $read1 $read2 \
        --outFileNamePrefix ${unmapped_dir}/${sample}/${sample}"." \
        --readFilesCommand zcat \
        #--limitGenomeGenerateRAM 358247314100 \
        --alignEndsType EndToEnd \
        --genomeDir $star_index \
        --alignIntronMax 1000000 \
        --outSAMunmapped Within \
        --runThreadN $threads \
        --twopassMode None \
        --outSAMtype BAM SortedByCoordinate
  else
      echo "[ ] STAR alignment on the reference is already done"
  fi
  
  # Subset the BAM file output of STAR to take unmapped pairs of reads: the ones were both mates are unmapped
  if [[ ! -f "${unmapped_dir}/$sample/${sample}.Unmapped.out.bam" ]]; then  
      echo -e "\n[ ] Extracting non mapping entries ..."  
      samtools view -b -@ $threads -f 12 "${unmapped_dir}/$sample/${sample}.Aligned.sortedByCoord.out.bam" > "${unmapped_dir}/$sample/${sample}.Unmapped.out.bam" && echo -e "[✓] Done\n"
  fi
  
  
  # Sort the unmapped reads bam file w.r.t. read names, this must be done before the netx step
  if [[ ! -f "${unmapped_dir}/$sample/${sample}.Unmapped.out.nameSorted.bam" ]]; then
      echo "[ ] Sorting unmapped bam"
      samtools sort -@ $threads -n "${unmapped_dir}/$sample/${sample}.Unmapped.out.bam" -o "${unmapped_dir}/$sample/${sample}.Unmapped.out.nameSorted.bam" && echo -e "[✓] Done\n"
  fi
  
  
  # Obtain the fastq files from the BAM file
  echo "[ ] Convertin BAM entries in new fastq files"
  samtools fastq  "${unmapped_dir}/$sample/${sample}.Unmapped.out.nameSorted.bam" \
      -1 "${unmapped_dir}/$sample/${sample}_R1.Unmapped.out" \
      -2 "${unmapped_dir}/$sample/${sample}_R2.Unmapped.out" \
      -0 /dev/null -s /dev/null \
      -n \
      --threads $threads && echo -e "[✓] Done\n"
      
  echo "[ ] Zipping unmapped fastq files"
  gzip "${unmapped_dir}/$sample/${sample}_R1.Unmapped.out"
  gzip "${unmapped_dir}/$sample/${sample}_R2.Unmapped.out" && echo -e "[✓] Done\n"

  # remove intermediate files to avoind occupying too much space
  echo "[ ] Removing intermediate BAM files to avoid too much storage usage"
  rm -f "${unmapped_dir}/$sample/${sample}.Unmapped.out.bam"
  rm -f "${unmapped_dir}/$sample/${sample}.Unmapped.out.nameSorted.bam"
  rm -f "${unmapped_dir}/$sample/${sample}.Aligned.sortedByCoord.out.bam" && echo -e "[✓] Done\n"
  
  
  # Report time taken for current sample
  end=$(date +%s)
  duration=$((end - start))
  printf -v duration_hms '%02d:%02d:%02d' $((duration/3600)) $(((duration%3600)/60)) $((duration%60))
  
  echo -e "Total time taken for retrieving unmapped reads in $sample: $duration_hms\n\n"
  
done

echo -e "[✓] STEP 1: RETRIVEAL OF UNMAPPED READS IS DONE"
echo "################################################"
echo


###############################################################################################################################
# sTEP 2 : METAPHLAN RUN TO MAP READS ON METAPHLAN DATABASE (from /data/adalla/scripts/run_mpa_bowtie2-paired.sh)

echo
echo "#######################################"
echo "STEP 2: ALIGNMENT TO METAPHLAN DATABASE"

for mate1_gz in $( ls "$unmapped_dir"/*/*_R1.Unmapped.out.gz ); do

    # Get corresponding mate2.gz file
    mate2_gz="${mate1_gz/_R1/_R2}"
    
    # Extract sample name from path
    sample=$(basename "$(dirname "$mate1_gz")")
    
    out="$outdir/$sample"
    mkdir -p "$out"
    
    [[ -f "$out/${sample}.bam" || -f "$out/${sample}-sorted.bam" ]] && echo -e "[✓] Skipping $sample" && continue
    
    echo "----- Precessing: $sample -----"
      
    # Define temporary uncompressed files
    mate1="${mate1_gz%.gz}"  # Remove .gz extension
    mate2="${mate2_gz%.gz}"  # Remove .gz extension

    # Define output file names
    mapout="$out/${sample}.mapout.bz2"
    samout="$out/${sample}.sam"
    vscout="$out/profiled-${sample}.vscout"
    output_file="$out/profiled-${sample}.txt"

    #echo "Writing in $BOWTIE2OUT"
    #echo "Writing in $OUTPUT_FILE"
    
    # Unzip files
    echo "Extracting $mate1"
    gunzip -c "$mate1_gz" > "$mate1"
    echo "Extracting $mate2"
    gunzip -c "$mate2_gz" > "$mate2"

    echo "[ ] Running MetaPhlAn"
    
    # Paired run of MetaPhlAn
    metaphlan ${mate1},${mate2} \
        --input_type fastq \
        --mapout "$mapout" \
        --nproc $threads \
        --profile_vsc \
        --vsc_out $vscout \
        --db_dir "$mpa_db" \
        --index "$mpa_index" \
        --samout "$samout" \
        --stat_q 0.01 \
        --force \
        --read_min_len 50 \
        --stat tavg_l \
        -t rel_ab_w_read_stats \
        -o "$output_file"

    # Remove uncompressed files after processing
    rm -f "$mate1"
    rm -f "$mate2"
    rm -f "$mapout"
    
    bam_file="${out}/${sample}.bam"
    samtools view -bS --threads $threads -o "$bam_file" "$samout"
    rm -f "$samout"
    
    echo "[✓] Processed $sample"
    echo
done

echo "[✓] All samples processed!"

echo "[✓] STEP 2: ALIGNMENT TO METAPHLAN DATABASE IS DONE"
echo "###################################################"
echo


########################################################################################################################
# STEP 3 : BLAST MAPPED READS ON HOMO SAPIENS AND FILTER THE ONES MAPPING (from /data/adalla/scripts/blast-onSapiens.sh)


echo
echo "##################################################"
echo "STEP 3: FILTERING OF HOMO SAPIENS READS WITH BLAST"

# First thing to do is retrieving the unique_species_file.txt that will be used in downstream analysis
python scripts/parse-metaphlan-out.py $wd $outdir "" 0

# Retreive all the markers associated to each species found in the dataset
mkdir -p $wd/species-markers

marker_info="$mpa_db/${mpa_index}_marker_info.txt"

if  [[ $(ls species-markers | wc -l) -ne $(wc -l < "$wd/preBLAST_species_file.tsv") ]]; then 
    echo "[*] Extracting markers for the found species"
    
    cut -f1 "$wd/preBLAST_species_file.tsv" | while read -r species; do
        # Output species markers file path
        species_markers_file="$wd/species-markers/markers-${species}.txt"
        
        # Grep and write to individual file
        grep "$species" "$marker_info" | cut -d$'\t' -f1 > "$species_markers_file"
    done
else
    echo -e "[✓] Markers file for each species already present\n"
fi



for sample_dir in "$outdir"/*; do

    # Start counting time
    start=$(date +%s)
    
    # Skip if not a directory
    [[ -d "$sample_dir" ]] || continue

    sample=$(basename "$sample_dir")
    
    samout="${sample_dir}/${sample}.sam"
    bam_file="${sample_dir}/${sample}.bam"
    sorted_bam="${sample_dir}/${sample}-sorted.bam"
    mkdir -p "${sample_dir}/BLAST-out"
    mapped_sorted_bam="${sample_dir}/BLAST-out/${sample}_mapped_reads.bam"
    mapped_filtered_sorted_bam="${sample_dir}/BLAST-out/${sample}_mapped_single_reads.bam"
    mapped_filtered_sorted_fasta="${sample_dir}/BLAST-out/${sample}_mapped_single_reads.fasta"
    blastout="${sample_dir}/BLAST-out/BLAST_${sample}_onSapiens.tsv"
    blast_best_hits="${sample_dir}/BLAST-out/BLAST_${sample}_best_hits.tsv"
    
    
    # [[ -f "$mpa_input_bam" ]] && echo -e "[✓] Skipping $sample" && continue

    # Since the pipeline support the creation of different mpa_input_bam with different s.i. 
    # i firstly check which ones are missing and which are already computed
    missing=()
    for si in "${si_thr[@]}"; do
        [[ ! -f "${sample_dir}/BLAST-out/${sample}_mpa-input_${si}.bam" ]] && missing+=("$si")
    done
    
    if [[ ${#missing[@]} -eq 0 ]]; then
        echo "--- Skipping $sample ---"
        continue
    else
        echo "----- Processing: $sample -----"
    fi
    
    # --- Ensure sorted & indexed BAM exists ---
    if [[ -f "${sorted_bam}" ]]; then
        echo "[✓] Sorted BAM already present"
    else
        # Need to generate sorted BAM and index, so ensure we have BAM
        if [[ ! -f "$bam_file" || ! -s "$bam_file" ]]; then
            echo "[ ] Need to compute the BAM file from the SAM"
            if [[ -f "$samout" ]]; then
                echo "[*] Converting SAM to BAM..."
                samtools view -bS --threads $threads -o "$bam_file" "$samout"
                echo "[*] Sorting BAM..."
                samtools sort --threads $threads -o "$sorted_bam" "$bam_file"
            else
                echo "[!] SAM file missing: $samout"
                echo "[!] MetaPhlAn run in the previous step did not end properly"
                continue  # Can't proceed without BAM
            fi
        else
            echo "[✓] BAM exists: $BAM_FILE"
            echo "[*] Sorting BAM..."
            samtools sort --threads $threads -o "$sorted_bam" "$bam_file"
        fi
    fi

    echo "[*] Indexing sorted BAM..."
    samtools index "$sorted_bam"
    

    if [[ ! -f "$mapped_sorted_bam" ]]; then
        echo "[*] Extracting reads mapping to some microbial sequences..."
        samtools view -h -b -F 2308 --threads $threads "$sorted_bam" > "$mapped_sorted_bam"
    
            # this flag is taken from https://www.biostars.org/p/138116/
            # Excludes: 
            # - 4 = unmapped reads
            # - 256 = secondary alignments = I'm not interested in them, if a reads aligns it is enough for me to blast it on human
            # - 2048 = supplementary alignments = I'm not interested in them, if a reads aligns it is enough for me to blast it on human
            # Actually in my BAM file there aren't unmapped reads because i used --no-unal option of bowtie2 
    else
        echo "[✓] $(basename $mapped_sorted_bam) already present"
    fi
    
    # To reduce the computational time of running blast, I do not take both reads of a pair (if they are both aligned), 
    # but i subset them so that i'm taking only single ends basically. This is because the reads comes from the same 
    # transcript, so i expect that the 2 mates would blast on the same sequence. Therefore having only one of them is
    # enough for me to blast
    if [[ ! -f "$mapped_filtered_sorted_bam" ]]; then
        echo "[*] Extracting single reads BAM file"
        samtools view -h "$mapped_sorted_bam" | { grep '^@'; samtools view "$mapped_sorted_bam" | sort -k1,1 -u; } | samtools view -b -o "$mapped_filtered_sorted_bam" -
            # This is the correct code because: 
              # samtools view -c DB_mapped_reads.bam  ->  120144
              # samtools view DB_mapped_reads.bam | cut -f1 | sort -u | wc -l  ->  103130
              # samtools view -c DB_mapped_single_reads.bam  ->  103130
    else
        echo "[✓] $(basename $mapped_filtered_sorted_bam) already present"
    fi
    
    tmp1=$( samtools view -c $mapped_sorted_bam )
    tmp2=$( samtools view -c $mapped_filtered_sorted_bam )
    diff=$(awk "BEGIN {printf \"%.2f\", 100 - ($tmp2 / $tmp1 * 100)}")
    echo "${sample}_mapped_reads.bam = $tmp1 reads"
    echo "${sample}_mapped_single_reads.bam = $tmp2"
    echo "About ${diff}% of the reads were saved"
    # echo "${sample},${diff}%" >> "$wd/saved-reads.csv"
    
    # If there is no output FASTA file of the reads to be blasted
    if [[ ! -f "$mapped_filtered_sorted_fasta"  || ! -s "$mapped_filtered_sorted_fasta" ]]; then
        echo "[*] Converting BAM to FASTA..."
        samtools fasta --threads $threads "$mapped_filtered_sorted_bam" > "$mapped_filtered_sorted_fasta"
    else
        echo "[✓] $(basename $mapped_filtered_sorted_fasta) already present"
    fi

    
    # Now i run BLAST looking for only matches on homo-sapiens. If i find them then I will remove that read in case of sequence identity of 100%
    echo "[*] Running megablast ..."
    blastn -query $mapped_filtered_sorted_fasta -db $blast_dir/$blast_index -task megablast -taxids 9606 -num_threads $threads -outfmt '6 qseqid sseqid qstart quend sstart send evalue bitscore length pident nident mismatch gapopen staxid ssciname sstrand' > "$blastout"

    [[ -f $mapped_filtered_sorted_fasta ]] && rm $mapped_filtered_sorted_fasta
    [[ -f $mapped_filtered_sorted_bam ]] && rm $mapped_filtered_sorted_bam
    [[ -f $mapped_sorted_bam ]] && rm $mapped_sorted_bam
    [[ -f $samout ]] && rm $samout
    [[ -f $bam_file ]] && rm $bam_file
    [[ -f "${sorted_bam}.bai" ]] && rm "${sorted_bam}.bai"
    

    end=$(date +%s)
    duration=$((end - start))
    printf -v duration_hms '%02d:%02d:%02d' $((duration/3600)) $(((duration%3600)/60)) $((duration%60))
    
    echo -e "Time taken for blasting $sample: $duration_hms\n"
    
    # Now we have to remove from the original bam file the reads that blasted on homo sapiens, and run the 
    # abundance computation step of Metaphlan on them, giving the sam file as input
    
    # We allow multiple thresholds, producing an output bam file each
    echo "[*] Removing reads blasted on homo sapiens"
    
    # Since it would be nice to store the $blastout file so that we don't have to 
    # reblast every sample anytime a new s.i.threshold is given in input, but they 
    # may be too big, we store only the top results w.r.t. s.i. for each read
    sort -k1,1 -k9,9nr "$blastout" | awk '!seen[$1]++' > "$blast_best_hits"
    
    for si in "${missing[@]}"; do
        echo -e "\n[*]At ${si}% of sequence identity ..."
        mpa_input_bam="${sample_dir}/BLAST-out/${sample}_mpa-input_${si}.bam"
        samtools view -h "$sorted_bam" | grep -v -F -f <(awk -v si_threshold="$si" '$9 >= si_threshold {sub(/\/.*/, "", $1); print $1}' "$blast_best_hits") | samtools view -Sb - > "$mpa_input_bam"
        
        # Legacy version with the blastout filtering and bam subsettin in one line:
        #samtools view -h "$sorted_bam" | grep -v -F -f <(sort -k1,1 -k9,9nr "$blastout" | awk -v si_threshold="$si" '!seen[$1]++ && $9 >= si_threshold {sub(/\/.*/, "", $1); print $1}') | samtools view -Sb - > $mpa_input_bam
    done

    rm -f $blastout
    
    end=$(date +%s)
    duration=$((end - start))
    printf -v duration_hms '%02d:%02d:%02d' $((duration/3600)) $(((duration%3600)/60)) $((duration%60))
    
    echo -e "Total time taken for STEP 3 on $sample: $duration_hms\n"
    echo "[✓] Processed $sample"
    echo
done

echo "[✓] All samples blasted!"    

    
echo "STEP 3: FILTERING OF HOMO SAPIENS READS WITH BLAST IS DONE"
echo "##########################################################"
echo


####################################################
# STEP 4 : FINAL METAPHLAN RUN ON THE FILTERED READS

echo
echo "##########################################################################"
echo "STEP 4: FINAL METAPHLAN RUN ON THE FILTERED READS AND COVERAGE COMPUTATION"


for sample_dir in "$outdir"/*; do

    # Skip if not a directory
    [[ -d "$sample_dir" ]] || continue
    
    # Completeness check
    sample=$(basename "$sample_dir")
    
    missing=()
    for si in "${si_thr[@]}"; do
        [[ ! -f "${sample_dir}/${sample}-final_${si}.txt" ]] && missing+=("$si")
    done
    
    if [[ ${#missing[@]} -eq 0 ]]; then
        echo "--- Skipping $sample ---"
        continue
    else
        echo -e "\n----- Processing: $sample -----"
    fi
    

    for si in "${missing[@]}"; do
        
        echo -e "[*]\tAt ${si}% of sequence identity ..."
        
        mpa_input_bam="${sample_dir}/BLAST-out/${sample}_mpa-input_${si}.bam"
        mpa_input_sam="${sample_dir}/BLAST-out/${sample}_mpa-input_${si}.sam"
        num_reads=$(  samtools view -c $mpa_input_bam )
    
        echo "reads in input to metaphlan: $num_reads"
        
        samtools view -h $mpa_input_bam > $mpa_input_sam
        
        metaphlan $mpa_input_sam \
            --input_type sam \
            --bowtie2db $mpa_db \
            --index $mpa_index \
            --stat_q 0.01 \
            --force \
            --stat tavg_l \
            -t rel_ab_w_read_stats \
            -o "${sample_dir}/${sample}-final_${si}.txt" \
            --nreads $num_reads
        
        [[ -f $mpa_input_sam ]] && rm $mpa_input_sam
        
        
        mkdir -p "$sample_dir/new-coverage"
        coverage_file="$sample_dir/new-coverage/${sample}-coverage-postBLAST_${si}.tsv"
        
        echo "[*] Computing coverage for $sample at s.i. ${si}%"
        samtools coverage -o "$coverage_file" $mpa_input_bam
            
    done

    echo "[✓] Finished processing $sample"
    echo
    
done

echo "[✓] All samples processed!"

echo "[✓] STEP 4: FINAL METAPHLAN RUN ON THE FILTERED READS AND COVERAGE COMPUTATION IS DONE"
echo "###############################à#######################################################"
echo


################################################
# STEP 5 : PLOTTING RESULTS AND COVERAGE STUDIES

echo "##############################################"
echo "STEP 5 : PLOTTING RESULTS AND COVERAGE STUDIES"


for si in "${si_thr[@]}"; do
    
    echo -e "\n\n----------------------------------------"
    echo -e "Results for sequence identity threshold at $si%"
    
    # Retrieve the final_species_file.tsv to compute the dataset report
    # Input: MetaPhlAn output files
    # Output: postBLAST_species_file_${si}.tsv
    python scripts/parse-metaphlan-out.py $wd $outdir $si 1

    
    for cov in "${coverage[@]}"; do
        
        echo -e "\n----- Results for sequence identity threshold at ${si}% with coverage ${cov}"
    
        # Apply the coverage threshold on the final_species_file.tsv
        # Input: postBLAST_species_file_${si}.tsv
        # Output: filtered_inspections_file_${si}_${coverage}.tsv
        python scripts/threshold-sample-report.py $wd $outdir $si $cov
        
        
        # Plot the results
        python scripts/new-plots.py $wd $si $cov
        
        
        # Statistical analysis and other plots
        # Rscript scripts/statistics.R $wd $si $cov
        
    done
done

# Plot the result to choose the threshold for the sample-report
# python scripts/choose-threshold.py "$wd" "$outdir" 

# Plot the sample-reports
# python scripts/sample-report.py "$wd" "$si_thr" "$plots_dir"


echo
echo "*************************************************************************"
echo "********************* PIPELINE ENDED SUCCESSFULLY ***********************"
echo "*************************************************************************"
echo
echo

    
    
