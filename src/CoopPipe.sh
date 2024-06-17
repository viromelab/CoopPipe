#!/bin/bash
#
RUN="0";
INSTALL="0";
INSTALL_SPECIFIC="0";
INSTALL_CONDA="0";
THREADS="4";
MEMORY="28";
#
CORONASPADES="0";
HAPLOFLOW="0";
LAZYPIPE="0";
METASPADES="0";
METAVIRALSPADES="0";
PEHAPLO="0";
QURE="0";
QVG="0";
SPADES="0";
SSAKE="0";
TRACESPIPE="0";
TRACESPIPELITE="0";
VIRGENA="0";
VISPA="0";
VPIPE="0";
#
VIRGENA_TIMEOUT="15";
#
FALCON_META="1";
KRAKEN2="0";
CENTRIFUGE="0";
#
BWA="0";
BOWTIE2="1";
#
MAFFT="1";
MUSCLE="0";
#
EMBOSS="0";
SELF="1";
#
READS1="";
READS2="";
CURR_PATH="$(pwd)";
TOOL_PATH="$(pwd)/HVRS/src";
OUTPUT="$(pwd)/out_analysis";
DATABASE="VDB.mfa";
declare -a VIRUSES_AVAILABLE=("B19V" "BuV" "CuV" "HBoV" "AAV" "BKPyV" "JCPyV" "KIPyV"
                    "WUPyV" "MCPyV" "HPyV6" "HPyV7" "TSPyV" "HPyV9" "MWPyV"
                    "STLPyV" "HPyV12" "NJPyV" "LIPyV" "SV40" "TTV" "TTVmid"
                    "TTVmin" "HAV" "HBV" "HCV" "HDV" "HEV" "SENV" "HPV2"
                    "HPV6" "HPV11" "HPV16" "HPV18" "HPV31" "HPV39" "HPV45"
                    "HPV51" "HPV56" "HPV58" "HPV59" "HPV68" "HPV77" "HSV-1"
                    "HSV-2" "VZV" "EBV" "HCMV" "HHV6" "HHV7" "KSHV" "ReDoV"
                    "VARV" "MPXV" "EV" "SARS2" "HERV" "MT");
declare -a VIRUSES;
FORCE_REFERENCES="0";
REFERENCES_DIR="";
AWK="0";
#
RESULT="0";
#
TOP="150";
################################################################################
#
SHOW_MENU () {
  echo " ------------------------------------------------------------------ ";
  echo "                                                                    ";
  echo " CoopPipe.sh : CoopPipe version v0.1                                ";
  echo "                                                                    ";
  echo " WIP version of CoopPipe. CoopPipe is a reconstruction pipeline     ";
  echo " that is capable of reconstructing human viral genomes, identifying ";
  echo " the viruses present and evaluating the reconstruction done.        ";  
  echo "                                                                    ";
  echo " Program options -------------------------------------------------- ";
  echo "                                                                    ";
  echo " -h, --help                    Show this,                           ";
  echo " --miniconda                   Install Miniconda,                   ";
  echo " -i, --install                 Install all tools (w/ conda),        ";
  echo " -is, --install-specific       Install only the tools specific to   ";
  echo "                               CoopPipe (use exclusively if all     ";
  echo "                               reconstruction tools have been       ";
  echo "                               installed using HVRS),               ";
  echo "                                                                    ";
  echo " --coronaspades                Reconstruction using coronaSPAdes,   ";
  echo " --haploflow                   Reconstruction using Haploflow,      ";
  echo " --lazypipe                    Reconstruction using LAZYPIPE,       ";
  echo " --metaspades                  Reconstruction using metaSPAdes,     ";
  echo " --metaviralspades             Reconstruction using metaviralSPAdes,";
  echo " --pehaplo                     Reconstruction using PEHaplo,        ";
  echo " --qure                        Reconstruction using QuRe,           ";
  echo " --qvg                         Reconstruction using QVG,            ";
  echo " --spades                      Reconstruction using SPAdes,         ";
  echo " --ssake                       Reconstruction using SSAKE,          ";
  echo " --tracespipe                  Reconstruction using TRACESPipe,     ";
  echo " --tracespipelite              Reconstruction using TRACESPipeLite, ";
  echo " --virgena                     Reconstruction using VirGenA,        ";
  echo " --vispa                       Reconstruction using ViSpA,          ";
  echo " --vpipe                       Reconstruction using V-pipe,         ";
  echo "                                                                    ";
  echo " --all                         Reconstruction using all tools,      ";
  echo "                                                                    ";
  echo " --virgena-timeout             Maximum time used by VirGenA         "; 
  echo "                               to reconstruct with each reference,  ";
  echo "                                                                    ";
  echo " --top_falcon  <INT>           Maximum number of references retrived";
  echo "                               by FALCON-meta,                      ";
  echo "                                                                    ";
  echo " -t  <INT>, --threads <INT>    Number of threads,                   ";
  echo " -m  <INT>, --memory <INT>     Maximum of RAM available,            ";
  echo " --tools  <STR>                Path to dir where reconstruction     ";
  echo "                               tools are installed,                 ";
  echo " --classifier  <STR>           Classifier to be used.               ";
  echo "                               Options: falcon, kraken2, centrifuge ";
  echo "                                                                    ";
  echo " --align  <STR>                Aligner to be used.                  ";
  echo "                               Options: bowtie2, bwa                ";
  echo "                                                                    ";
  echo " --msa  <STR>                  MSA to be used.                      ";
  echo "                               Options: mafft, muscle               ";
  echo "                                                                    ";
  echo " --consensus  <STR>            Consensus script to be used.         ";
  echo "                               Options: self, emboss                ";
  echo "                                                                    ";
  echo " --awk  <STR>                  Generate consensus with AWK.         ";
  echo " --emboss  <STR>               Generate consensus with EMBOSS.      ";
  echo "                                                                    ";
  echo "                                                                    ";
  echo " -o  <STR>, --output <STR>     Output folder name,                  ";
  echo "                                                                    ";
  echo " -r1 <STR>, --reads1 <STR>     FASTQ reads (forward),               ";
  echo " -r2 <STR>, --reads2 <STR>     FASTQ reads (reverse),               ";
  echo "                                                                    ";
  echo " Examples --------------------------------------------------------- ";
  echo "                                                                    ";
  echo " - Install Miniconda                                                ";
  echo "  ./CoopPipe.sh --miniconda                                         ";
  echo "                                                                    ";
  echo " - Install tools                                                    ";
  echo "  ./CoopPipe.sh -i --threads 4                                      ";
  echo "                                                                    "; 
  echo " - Running example                                                  ";
  echo "  ./CoopPipe.sh -r1 reads_forward.fq -r2 reads_reverse.fq  \\       ";
  echo "    --output out_analysis --threads 4 --memory 28          \\       ";
  echo "    --tools HVRS/src --classifier falcon --all                      ";
  echo "                                                                    ";
  echo " ------------------------------------------------------------------ ";
  }
#
################################################################################
#
CHECK_INPUT () {
  FILE=$1
  if [ -f "$FILE" ];
    then
    echo "Input filename: $FILE"
    else
    echo -e "\e[31mERROR: input file not found ($FILE)!\e[0m";
    SHOW_MENU;
    exit;
    fi
  }
#
################################################################################
#
PROGRAM_EXISTS () {
  printf "Checking $1 ... ";
  if ! [ -x "$(command -v $1)" ];
    then
    echo -e "\e[41mERROR\e[49m: $1 is not installed." >&2;
    echo -e "\e[42mTIP\e[49m: Try: ./CoopPipe.sh --install" >&2;
    exit 1;
    else
    echo -e "\e[42mSUCCESS!\e[49m";
    fi
  }
#
################################################################################
#
CHECK_PROGRAMS () {
  PROGRAM_EXISTS "AdapterRemoval";
  PROGRAM_EXISTS "gto_fasta_extract_read_by_pattern";
  PROGRAM_EXISTS "gto_filter";
  PROGRAM_EXISTS "FALCON";
  PROGRAM_EXISTS "bwa";
  PROGRAM_EXISTS "samtools";
  PROGRAM_EXISTS "bcftools";
  PROGRAM_EXISTS "bedops";
  PROGRAM_EXISTS "bedtools";
  PROGRAM_EXISTS "igv";
  PROGRAM_EXISTS "tabix";
  }
#
################################################################################
#
#creates a fasta file with paired reads
CREATE_PAIRED_FA_FILES () { 
  printf "Creating fasta files from .fq files\n\n"	    
  sed -n '1~4s/^@/>/p;2~4p' $READS1 > tmp_new_1.fa
  sed -n '1~4s/^@/>/p;2~4p' $READS2 > tmp_new_2.fa
  cat tmp_new_*.fa > input.fasta
  perl -pe 's/[\r\n]+/;/g; s/>/\n>/g' input.fasta | sort -t"[" -k2,2V | sed 's/;/\n/g' | sed '/^$/d' > tmp.txt
    rm tmp.txt
  mv input.fasta paired.fa
}

#
## RUN VIRAL METAGENOMIC COMPOSITION =========================================
#
CLASSIFY_INPUT () {
  OUTPUT_DIR=$1
  INPUT_FILE=$2
    
  aux_OUTPUT_DIR="$(echo $OUTPUT_DIR | awk -F/ '{print $NF}')" 
  lzma -d -k -f VDB.mfa.lzma 
  
  if [[ "$FALCON_META" -eq "1" ]];
    then
    printf "Classifying with FALCON-meta.\n\n"
    conda activate falcon
  
    #lzma -d $DATABASE.lzma
    lzma -d VDB.fa.lzma
  
    FALCON -v -n $THREADS -t $TOP -F -m 6:1:1:0/0 -m 13:50:1:0/0 -m 19:500:1:5/10 -g 0.85 -c 10 -x top-metagenomics.csv $INPUT_FILE $DATABASE
    #
    ## GET HIGHEST SIMILAR REFERENCE =============================================
    #
    rm -f best-viral-metagenomics.txt
    for VIRUS in "${VIRUSES_AVAILABLE[@]}"
      do
      printf "%s\t" "$VIRUS" >> best-viral-metagenomics.txt;
      #
      RESULT=`cat top-metagenomics.csv | grep -a -f IDS/ID-$VIRUS.ids \
      | awk '{ if($3 > 0 && $2 > 200 && $2 < 9000000) print $3"\t"$4; }' \
      | head -n 1 \
      | awk '{ print $1"\t"$2;}' \
      | sed "s/NC\_/NC-/" \
      | tr '_' '\t' \
      | awk '{ print $1"\t"$2;}'`;
      if [ -z "$RESULT" ]
        then
        echo -e "-\t-" >> best-viral-metagenomics.txt
        else
        echo "$RESULT" | sed "s/NC-/NC\_/" >> best-viral-metagenomics.txt
        fi
      done
  
    rm -rf $OUTPUT_DIR
    mkdir $OUTPUT_DIR
    while IFS= read -r line
    do
      PERC="$(cut -d'	' -f2 <<< $line)"
      if [ "$PERC" != "-" ];
        then
        NAME="$(cut -d'	' -f1 <<< $line)"
        ID="$(cut -d'	' -f3 <<< $line)"
        gto_fasta_extract_read_by_pattern -p $ID < VDB.fa | sed 's|/||g' > $NAME.fa
        mv $NAME.fa $OUTPUT_DIR
      fi 
    done < best-viral-metagenomics.txt
   
    mv best-viral-metagenomics.txt $OUTPUT_DIR
    mv top-metagenomics.csv $OUTPUT_DIR
    conda activate base
  
  
  elif [[ "$KRAKEN2" -eq "1" ]];
    then 
    printf "Classifying with Kraken2.\n\n"
    eval "$(conda shell.bash hook)"
    rm -rf $OUTPUT_DIR 
    mkdir $OUTPUT_DIR

    conda activate kraken
    kraken2 --threads $THREADS --db krakendb/ --output output.txt --report report.txt $INPUT_FILE
    conda activate base
    echo "" > list_species_refs.txt
    while IFS= read -r line
      do
      RANK="$(cut -d'	' -f4 <<< $line)"      

      if [[ "$RANK" == "S" || "$RANK" == "S1" || "$RANK" == "S2" ]];
        then        
        SPECIES_NAME="$(cut -d'	' -f6 <<< $line)"
        SPECIES_NAME=$(echo "$SPECIES_NAME" | tr -s " ")
        LIST_POSSIBILITIES=$(cat VDB.fa | grep -i "$SPECIES_NAME" )
        #echo "$line" >> list_species_refs.txt       
        for VIRUS in "${VIRUSES_AVAILABLE[@]}"
          do          
          RESULT=$(echo "$LIST_POSSIBILITIES" | grep -i -a -f IDS/ID-$VIRUS.ids | head -n 1)
          #echo "$RESULT" >> list_species_refs.txt 
          AUX="$(cut -d' ' -f1 <<< $RESULT)"
          ID="$(cut -d'>' -f2 <<< $AUX)"
          if [[ "$ID" ]];
            then
            printf "$ID \n"
            gto_fasta_extract_read_by_pattern -p $ID < VDB.fa | sed 's|/||g' > $VIRUS.fa
            mv $VIRUS.fa $OUTPUT_DIR
          fi
        done
      fi 
    done < report.txt
    mv report.txt $OUTPUT_DIR
    mv output.txt $OUTPUT_DIR
    sed -i '/^$/d' list_species_refs.txt     
    
  elif [[ "$CENTRIFUGE" -eq "1" ]];
    then
    printf "Classifying with centrifuge.\n\n"
    eval "$(conda shell.bash hook)"
    rm -rf $OUTPUT_DIR 
    mkdir $OUTPUT_DIR

    conda activate centrifuge
    centrifuge -f -x abv $INPUT_FILE
    conda activate base
    
    while IFS= read -r line
      do
      SPECIES_NAME="$(cut -d'	' -f1 <<< $line)"     
      LIST_POSSIBILITIES=$(cat VDB.fa | grep -i "$SPECIES_NAME" )
      #echo "$line" >> list_species_refs.txt       
      for VIRUS in "${VIRUSES_AVAILABLE[@]}"
        do          
        RESULT=$(echo "$LIST_POSSIBILITIES" | grep -i -a -f IDS/ID-$VIRUS.ids | head -n 1)
        #echo "$RESULT" >> list_species_refs.txt 
        AUX="$(cut -d' ' -f1 <<< $RESULT)"
        ID="$(cut -d'>' -f2 <<< $AUX)"
        if [[ "$ID" ]];
          then
          printf "$ID \n"
          gto_fasta_extract_read_by_pattern -p $ID < VDB.fa | sed 's|/||g' > $VIRUS.fa
          mv $VIRUS.fa $OUTPUT_DIR
        fi
      done   
    done < centrifuge_report.tsv
    mv centrifuge_report.tsv $OUTPUT_DIR
  fi
  #centrifuge-build -p 4 --conversion-table taxonomy/seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp VDB.fa abv

}

CLASSIFY_TRACES () {
  FILE=$1
  OUT_DIR=$2
  
  mkdir $OUT_DIR
  
  SAMPLES=( $(cat $FILE | tr "\t" "~" | cut -d"~" -f3 ) )

  for sample in "${SAMPLES[@]}"
      do
      if [[ "$sample" != "-" ]];
        then
        NAME=$( cat $FILE | grep $sample| tr "\t" "~" | cut -d"~" -f1  ) 
        if [[ "$NAME" == "B19" ]];
          then
          NAME="B19V"
        elif [[ "$NAME" == "HPV" ]];
          then
          NAME="HPV68"
        elif [[ "$NAME" == "HV3" ]];
          then
          NAME="VZV"
        elif [[ "$NAME" == "POLY5" ]];
          then
          NAME="MCPyV"
        fi
        gto_fasta_extract_read_by_pattern -p "$sample" < $CURR_PATH/VDB.fa > "$NAME.fa"
        mv "$NAME.fa" "$OUT_DIR"
      fi
  done
  mv $FILE $OUT_DIR

}
#
## RUN VIRAL METAGENOMIC COMPOSITION =========================================
#
ALIGN_TO_REF () {
  REF_DIR=$1
  INPUT=$2
  NAME_TOOL=$3
  #
  eval "$(conda shell.bash hook)"
  conda activate align
  #
  cd $REF_DIR
  #
  if [[ "$BWA" -eq "1" ]];
    then
    printf "Creating consensus with BWA.\n\n"
    for file in $(ls *.fa);
     do
      printf "$file \n\n"
      aux_file="$(cut -d'.' -f1 <<< $file)"

      cp $file SPECIFIC.fa

      bwa index SPECIFIC.fa
      bwa aln -t $THREADS SPECIFIC.fa $INPUT > SPECIFIC-READS.sai
      bwa samse SPECIFIC.fa SPECIFIC-READS.sai $INPUT > SPECIFIC-READS.sam
      conda activate sambcf
      samtools view -bSh SPECIFIC-READS.sam > SPECIFIC-READS.bam;
      samtools view -bh -F4 SPECIFIC-READS.bam > FIL-SPECIFIC-READS.bam;
      samtools sort -o SORT-FIL-SPECIFIC-READS.bam FIL-SPECIFIC-READS.bam;
      samtools rmdup -s SORT-FIL-SPECIFIC-READS.bam RD-SORT-FIL-SPECIFIC-READS.bam;
      samtools index -b RD-SORT-FIL-SPECIFIC-READS.bam RD-SORT-FIL-SPECIFIC-READS.bam.bai
      #
      bedtools genomecov -ibam RD-SORT-FIL-SPECIFIC-READS.bam -bga > SPECIFIC-coverage.bed
      awk '$4 < 1' SPECIFIC-coverage.bed > SPECIFIC-zero-coverage.bed
      samtools faidx SPECIFIC.fa
      bcftools mpileup -Ou -f SPECIFIC.fa RD-SORT-FIL-SPECIFIC-READS.bam \
      | bcftools call --ploidy 1 -P 9.9e-1 -mv -Oz -o SPECIFIC-calls.vcf.gz
      bcftools index SPECIFIC-calls.vcf.gz
      bcftools norm -f SPECIFIC.fa SPECIFIC-calls.vcf.gz -Oz -o SPECIFIC-calls.norm.vcf.gz
      bcftools filter --IndelGap 5 SPECIFIC-calls.norm.vcf.gz -Oz -o SPECIFIC-calls.norm.flt-indels.vcf.gz
      zcat SPECIFIC-calls.norm.flt-indels.vcf.gz | vcf2bed --snvs > SPECIFIC-calls.bed
      tabix -f SPECIFIC-calls.norm.flt-indels.vcf.gz
      bcftools consensus -m SPECIFIC-zero-coverage.bed -f SPECIFIC.fa SPECIFIC-calls.norm.flt-indels.vcf.gz > SPECIFIC-consensus.fa
      tail -n +2 SPECIFIC-consensus.fa > SPECIFIC-TMP_FILE.xki
      echo ">SPECIFIC consensus (REF: $SPECIFIC) [$NAME_TOOL]" > SPECIFIC-consensus.fa
      cat SPECIFIC-TMP_FILE.xki >> SPECIFIC-consensus.fa
      mv SPECIFIC-consensus.fa $NAME_TOOL-$aux_file-consensus.fa 
      mv $NAME_TOOL-$aux_file-consensus.fa $OUTPUT/consensus
      conda activate align
      rm -f SPECIFIC-TMP_FILE.xki;
      rm *SPECIFIC*
    done
  elif [[ "$BOWTIE2" -eq "1" ]];
    then
    printf "Creating consensus with Bowtie2.\n\n"
    for file in $(ls *.fa);
     do
      printf "$file \n\n"
      aux_file="$(cut -d'.' -f1 <<< $file)"
      

      cp $file SPECIFIC.fa
      bowtie2-build SPECIFIC.fa host_DB
      bowtie2 -p 1 -x host_DB -f $INPUT > SPECIFIC-READS.sam #2>> report_stderr.txt
      conda activate sambcf
      samtools view -bSh SPECIFIC-READS.sam > SPECIFIC-READS.bam;
      samtools view -bh -F4 SPECIFIC-READS.bam > FIL-SPECIFIC-READS.bam;
      samtools sort -o SORT-FIL-SPECIFIC-READS.bam FIL-SPECIFIC-READS.bam;
      samtools rmdup -s SORT-FIL-SPECIFIC-READS.bam RD-SORT-FIL-SPECIFIC-READS.bam;
      samtools index -b RD-SORT-FIL-SPECIFIC-READS.bam RD-SORT-FIL-SPECIFIC-READS.bam.bai
      #
      bedtools genomecov -ibam RD-SORT-FIL-SPECIFIC-READS.bam -bga > SPECIFIC-coverage.bed
      awk '$4 < 1' SPECIFIC-coverage.bed > SPECIFIC-zero-coverage.bed
      samtools faidx SPECIFIC.fa
      bcftools mpileup -Ou -f SPECIFIC.fa RD-SORT-FIL-SPECIFIC-READS.bam \
      | bcftools call --ploidy 1 -P 9.9e-1 -mv -Oz -o SPECIFIC-calls.vcf.gz
      bcftools index SPECIFIC-calls.vcf.gz
      bcftools norm -f SPECIFIC.fa SPECIFIC-calls.vcf.gz -Oz -o SPECIFIC-calls.norm.vcf.gz
      bcftools filter --IndelGap 5 SPECIFIC-calls.norm.vcf.gz -Oz -o SPECIFIC-calls.norm.flt-indels.vcf.gz
      zcat SPECIFIC-calls.norm.flt-indels.vcf.gz | vcf2bed --snvs > SPECIFIC-calls.bed
      tabix -f SPECIFIC-calls.norm.flt-indels.vcf.gz
      bcftools consensus -m SPECIFIC-zero-coverage.bed -f SPECIFIC.fa SPECIFIC-calls.norm.flt-indels.vcf.gz > SPECIFIC-consensus.fa
      tail -n +2 SPECIFIC-consensus.fa > SPECIFIC-TMP_FILE.xki
      echo ">SPECIFIC consensus (REF: $SPECIFIC) [$NAME_TOOL]" > SPECIFIC-consensus.fa
      cat SPECIFIC-TMP_FILE.xki >> SPECIFIC-consensus.fa
      mv SPECIFIC-consensus.fa $NAME_TOOL-$aux_file-consensus.fa 
      mv $NAME_TOOL-$aux_file-consensus.fa $OUTPUT/consensus
      conda activate align
      rm -f SPECIFIC-TMP_FILE.xki;
      rm *SPECIFIC*
      rm host_DB*
      done    
  fi

  conda activate base
  cd $TOOL_PATH


}

PERFORM_MULTIPLE_ALIGNMENT () {
  DIRECTORY=$1
  #
  eval "$(conda shell.bash hook)"
  #
  if [[ "$MAFFT" -eq "1" ]];
    then
    printf "Performing multiple alignment using MAFFT.\n\n"
    conda activate mafft
    cd $DIRECTORY
    for virus in "${VIRUSES[@]}"
      do
      aux_virus="$(cut -d'.' -f1 <<< $virus)"
      aux_virus="$(echo $aux_virus | awk -F/ '{print $NF}')"
      #ls
      cat *-$aux_virus-consensus.fa > $aux_virus-combined.fa
      mafft --auto $aux_virus-combined.fa > multifasta-$aux_virus.fa
      
    done
    cd $CURR_PATH
    conda activate base
  elif [[ "$MUSCLE" -eq "1" ]];
    then
    printf "Performing multiple alignment using MUSCLE.\n\n"
    conda activate muscle
    cd $DIRECTORY
    for virus in "${VIRUSES[@]}"
      do
      aux_virus="$(cut -d'.' -f1 <<< $virus)"
      aux_virus="$(echo $aux_virus | awk -F/ '{print $NF}')"
      #ls
      cat *-$aux_virus-consensus.fa > $aux_virus-combined.fa
      muscle -align $aux_virus-combined.fa -output multifasta-$aux_virus.fa
    done
    cd $CURR_PATH
    conda activate base
  fi
}

CREATE_FINAL_CONSENSUS () {
  DIRECTORY=$1
  #
  cd $DIRECTORY
  for file in $(ls multifasta*.fa);
    do
    eval "$(conda shell.bash hook)"
    tmp="$(cut -d'.' -f1 <<< $file)"
    name_vir="$(cut -d'-' -f2 <<< $tmp)"    
    #cons -sequence $file -outseq cooppipeemboss-$name_vir-consensus.fa
    #python3 $CURR_PATH/generate_consensus.py -i $file -o cooppipealt-$name_vir-consensus.fa
    if [[ "$EMBOSS" -eq "1" ]];
      then
      
      conda activate emboss
      cons -sequence $file -outseq cooppipe-$name_vir-consensus.fa
      conda activate base
    fi
    
    if [[ "$SELF" -eq "1" ]];
      then
      python3 $CURR_PATH/weighted_generate_consensus.py -i $file -v $name_vir -k 5 15 30 100 200 400 500
      python3 $CURR_PATH/generate_consensus.py -i new.fa -o cooppipe-$name_vir-consensus.fa
    fi
  done
  #mkdir combined
  #mv *-combined.fa combined
  #mkdir multifasta
  #mv multifasta-*.fa multifasta

  cd $CURR_PATH
  
}
#
check_installation () { 
  NAME_TOOL=$1
  RESULT=0
  $TOOL_PATH/Verification.sh --$NAME_TOOL > verif.txt
  if [[ $(wc -l < verif.txt ) -eq "1" ]] 
    then
    RESULT=1
  else
    RESULT=0
  fi
  #printf "$NAME_TOOL      $RESULT\n"
}
#
################################################################################
#
if [[ "$#" -lt 1 ]];
  then
  HELP=1;
  fi
#
POSITIONAL=();
#
while [[ $# -gt 0 ]]
  do
  i="$1";
  case $i in
    -h|--help|?)
      HELP=1;
      shift
    ;;
    --miniconda)
      INSTALL_CONDA=1;
      shift
    ;;
    -i|--install)
      INSTALL=1;
      INSTALL_SPECIFIC=1;
      shift
    ;;
    -is|--install-specific)
      INSTALL_SPECIFIC=1;
      shift
    ;;
      --coronaspades)
      check_installation coronaspades;
      CORONASPADES=$RESULT;
      shift
    ;;
    --haploflow)
      check_installation haploflow;
      HAPLOFLOW=$RESULT;
      shift
    ;;
    --lazypipe)
      check_installation lazypipe;
      LAZYPIPE=$RESULT;
      shift
    ;;
    --metaspades)
      check_installation metaspades;
      METASPADES=$RESULT;
      shift
    ;;
    --metaviralspades)
      check_installation metaviralspades;
      METAVIRALSPADES=$RESULT;
      shift
    ;;
    --pehaplo)
      check_installation pehaplo;
      PEHAPLO=$RESULT;
      shift
    ;;
    --qure)
      check_installation qure;
      QURE=$RESULT;
      shift
    ;;
    --qvg)
      check_installation qvg;
      QVG=$RESULT;
      shift
    ;;
    --spades)      
      check_installation spades;
      SPADES=$RESULT;
      printf "SPADES -> $RUN_SPADES\n\n"
      shift
    ;;
    --ssake)
      check_installation ssake;
      SSAKE=$RESULT;
      shift
    ;;
    --tracespipe)
      check_installation tracespipe;
      TRACESPIPE=$RESULT;
      shift
    ;;
    --tracespipelite)
      check_installation tracespipelite;
      TRACESPIPELITE=$RESULT;
      shift
    ;;
    --virgena)
      check_installation virgena;
      VIRGENA=$RESULT;
      shift
    ;;
    --vispa)
      check_installation vispa;
      VISPA=$RESULT;
      shift
    ;;
    --vpipe)
      check_installation vpipe;
      VPIPE=$RESULT;
      shift
    ;;
    --all)
      check_installation coronaspades;
      CORONASPADES=$RESULT;
      check_installation haploflow;
      HAPLOFLOW=$RESULT;
      check_installation lazypipe;
      LAZYPIPE=$RESULT;
      check_installation metaspades;
      METASPADES=$RESULT;
      check_installation metaviralspades;
      METAVIRALSPADES=$RESULT;
      check_installation pehaplo;
      PEHAPLO=$RESULT;
      check_installation qure;
      QURE=$RESULT;
      check_installation qvg;
      QVG=$RESULT;
      check_installation spades;
      SPADES=$RESULT;
      check_installation ssake;
      SSAKE=$RESULT;
      check_installation tracespipe;
      TRACESPIPE=$RESULT;
      check_installation tracespipelite;
      TRACESPIPELITE=$RESULT;
      check_installation virgena;
      VIRGENA=$RESULT;
      check_installation vispa;
      VISPA=$RESULT;
      check_installation vpipe;
      VPIPE=$RESULT;
      shift
    ;;
    --virgena-timeout)
      VIRGENA_TIMEOUT="$2";
      shift 2;
    ;;
    -t|--threads)
      THREADS="$2";
      shift 2;
    ;;
    -m|--memory)
      MEMORY="$2";
      shift 2;
    ;;
    --tools)
      TOOL_PATH="$2";
      shift 2;
    ;;
    --classifier)
      TMP="$2";
      TMP=$(echo "$TMP" | tr '[:upper:]' '[:lower:]')

      if [[ "$TMP" == "kraken2" ]]; 
        then
        KRAKEN2="1"
        FALCON_META="0" 
      elif [[ "$TMP" == "centrifuge" ]];
        then
        CENTRIFUGE="1"
        FALCON_META="0" 
      elif [[ "$TMP" == "falcon" ]];
        then
        FALCON_META="1"  
      else
        FORCE_REFERENCES="1";
        REFERENCES_DIR="$2";
      fi     
      shift 2;
    ;;
    --top_falcon)
      TOP="$2";
      shift 2;
    ;;
    --align)
      TMP="$2";
      TMP=$(echo "$TMP" | tr '[:upper:]' '[:lower:]')

      if [[ "$TMP" == "BWA" ]];
        then
        BWA="1"
        BOWTIE2="0"  
      fi     
      shift 2;
    ;;
    --msa)
      TMP="$2";
      TMP=$(echo "$TMP" | tr '[:upper:]' '[:lower:]')

      if [[ "$TMP" == "muscle" ]];
        then
        MUSCLE="1"
        MAFFT="0"  
      fi     
      shift 2;
    ;;
    --consensus)
      TMP="$2";
      TMP=$(echo "$TMP" | tr '[:upper:]' '[:lower:]')

      if [[ "$TMP" == "emboss" ]];
        then
        EMBOSS="1"
        SELF="0"  
      fi     
      shift 2;
    ;;
    --awk)
      OUTPUT="$(pwd)/$2";
      AWK="1"
      shift 2;
    ;;
    --emboss)
      OUTPUT="$(pwd)/$2";
      EMBOSS="1"
      SELF="0"  
      shift 2;
    ;;
    -o|--output)
      OUTPUT="$(pwd)/$2";
      shift 2;
    ;;
    -r1|-r|--input1|--reads|--reads1)
      READS1="$2";
      RUN=1;
      shift 2;
    ;;
    -r2|--input2|--reads2)
      READS2="$2";
      RUN=1;
      shift 2;
    ;;
    -*) # unknown option with small
    echo "Invalid arg ($1)!";
    echo "For help, try: CoopPipe.sh -h"
    exit 1;
    ;;
  esac
  done
#
set -- "${POSITIONAL[@]}" # restore positional parameters
#
################################################################################
#
if [[ "$HELP" -eq "1" ]];
  then
  SHOW_MENU;
  exit;
  fi
#
################################################################################
#
if [[ "$INSTALL_CONDA" -eq "1" ]];
  then
  #
  rm -rf Miniconda3-latest-Linux-x86_64.sh
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  chmod +x Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  printf "Please close and reopen this terminal \n\n"
  read a
  #
  fi
#
################################################################################
#
if [[ "$INSTALL_SPECIFIC" -eq "1" ]];
  then
  #
  eval "$(conda shell.bash hook)"
  sudo apt install git -y
  #rm -rf HVRS
  if [[ "$INSTALL" -eq "1" ]];
    then
    git clone https://github.com/mirakaya/HVRS.git
    cd HVRS/src/
    chmod +x *.sh
    ./Installation.sh --all
    cd ../../
  fi
  #FALCON-meta
  conda create -n falcon -y
  conda activate falcon 
  conda install -c cobilab -y gto
  conda install -c cobilab falcon -y
  conda activate base
  #kraken2
  conda create -n kraken -y
  conda activate kraken
  conda install kraken2 -y
  lzma -d VDB.fa.lzma 
  kraken2-build --db krakendb --download-taxonomy --use-ftp
  kraken2-build --db krakendb --add-to-library VDB.fa
  kraken2-build --db krakendb --build --threads $THREADS
  conda install -c cobilab -y gto
  
  #kraken2 --threads $THREADS --db /opt/storage2/db/kraken2/standard --output ERR2513180.output.txt --report ERR2513180.report.txt --paired ERR2513180_1.fastq.gz ERR2513180_2.fastq.gz
  conda activate base
  #centrifuge
  conda create -n centrifuge -y
  conda activate centrifuge
  conda install -c bioconda centrifuge -y
  conda install -c cobilab -y gto
  centrifuge-download -o taxonomy taxonomy
  cp krakendb/seqid2taxid.map taxonomy
  centrifuge-build -p $THREADS --conversion-table taxonomy/seqid2taxid.map --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp VDB.fa abv
  conda activate base
  #bwa and bowtie2
  conda create -n align -y
  conda activate align
  conda install -c bioconda bwa bowtie2 -y
  #bowtie2-build VDB.fa host_DB
  conda activate base
  #sam and bcftools
  conda create -n sambcf -y
  conda activate sambcf
  conda install -c bioconda bcftools -y
  conda install -c bioconda samtools -y
  conda install -c bioconda bedtools -y
  conda activate base
  #mafft 
  conda create -n mafft -y 
  conda activate mafft
  conda install -c bioconda mafft -y
  conda activate base
  #muscle 
  conda create -n muscle -y 
  conda activate muscle
  conda install -c bioconda muscle -y
  conda activate base
  #emboss
  conda create -n emboss -y 
  conda activate emboss
  conda install -c bioconda emboss -y
  conda activate base
  #
  ./HVRS/src/Verification.sh --tools
  #
  fi
#
################################################################################
#
if [[ "$RUN" -eq "1" ]];
  then
  eval "$(conda shell.bash hook)"
  #
  CHECK_INPUT $READS1
  CHECK_INPUT $READS2
  #
  printf "Running CoopPipe with $THREADS threads.\n\n"
  cp $READS1 .
  cp $READS2 .
  aux_READS1="$(echo $READS1 | awk -F/ '{print $NF}')"
  aux_READS2="$(echo $READS2 | awk -F/ '{print $NF}')"
  #
  READS1="$(pwd)/$(echo $READS1 | awk -F/ '{print $NF}')"
  READS2="$(pwd)/$(echo $READS2 | awk -F/ '{print $NF}')"  
  #
  printf "$(pwd)\n\n"
  CREATE_PAIRED_FA_FILES
  if [[ "$FORCE_REFERENCES" -eq "1" ]];
    then
    printf "References were forced. Path is $REFERENCES_DIR\n\n"
    cd $REFERENCES_DIR
    VIRUSES=($(ls -1 -d "$PWD/"*.fa)) 
    cd $CURR_PATH
  else
    CLASSIFY_INPUT references paired.fa
    cd references
    VIRUSES=($(ls -1 -d "$PWD/"*.fa)) 
    cd ..
  fi
  #
  cd $TOOL_PATH
  #
  rm -rf $OUTPUT
  mkdir $OUTPUT
  mkdir $OUTPUT/consensus
  #
  ## RUN TOOLS =================================================================
  #
if [[ "$CORONASPADES" -eq "1" ]];
    then
    printf "Reconstructing with coronaSPAdes\n\n"
    eval "$(conda shell.bash hook)"
    conda activate spades
    rm -rf coronaspades_reconstruction
    mkdir coronaspades_reconstruction
    cd coronaspades_reconstruction

    mkdir input
    rm -rf output
    mkdir output	
    cp $READS1 $READS2 input
    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o coronaspades-time.txt coronaspades.py -o output -1 input/$aux_READS1 -2 input/$aux_READS2 -t $THREADS -m $MEMORY 
    mv coronaspades-time.txt $OUTPUT
    mv output/raw_scaffolds.fasta output/coronaspades-reconstructed.fa
    cp output/coronaspades-reconstructed.fa $OUTPUT
    cd ..
    conda activate base
    rm -rf coronaspades_reconstruction
    cd $CURR_PATH
    CLASSIFY_INPUT coronaSPAdes $OUTPUT/coronaspades-reconstructed.fa
    mv coronaSPAdes $OUTPUT
    ALIGN_TO_REF $OUTPUT/coronaSPAdes $OUTPUT/coronaspades-reconstructed.fa coronaspades
    cd $TOOL_PATH
  fi
  #
  if [[ "$HAPLOFLOW" -eq "1" ]];
    then
    printf "Reconstructing with Haploflow\n\n"
    eval "$(conda shell.bash hook)"
    conda activate haploflow
    rm -rf haploflow_data
    mkdir haploflow_data
    cd haploflow_data
    cp $READS1 $READS2 .
    cat $aux_READS1 $aux_READS2 > paired.fq
    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o haploflow-time.txt haploflow --read-file paired.fq --out test --log test/log    
    mv haploflow-time.txt $OUTPUT
    mv test/contigs.fa test/haploflow-reconstructed.fa
    cp test/haploflow-reconstructed.fa $OUTPUT
    rm -rf test
    cd ..  
    conda activate base 
    rm -rf haploflow_data
    cd $CURR_PATH
    CLASSIFY_INPUT Haploflow $OUTPUT/haploflow-reconstructed.fa
    mv Haploflow $OUTPUT
    cd $TOOL_PATH
  fi
  #
  if [[ "$LAZYPIPE" -eq "1" ]];
    then
    printf "Reconstructing with Lazypipe\n\n"
    eval "$(conda shell.bash hook)"
    conda activate blast
    conda activate --stack lazypipe
    cd lazypipe
    export TM="$CONDA_PREFIX/share/trimmomatic"
    export BLASTDB="$(pwd)/BLASTDB/"
    export data="$(pwd)/data/"
    cp $READS1 $READS2 .
    rm -rf results
    mkdir results
    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o lazypipe-time.txt perl lazypipe.pl -1 $aux_READS1 -2 $aux_READS2 --pipe all,rep -v -t $THREADS -res output

    mv output/*/contigs.fa results/lazypipe-reconstructed.fa
    rm -rf output
    mkdir output
    
    cp results/lazypipe-reconstructed.fa $OUTPUT
    mv lazypipe-time.txt $OUTPUT
    rm -rf *.fq
    rm -rf results*
    cd ..
    cd $CURR_PATH
    CLASSIFY_INPUT LAZYPIPE $OUTPUT/lazypipe-reconstructed.fa
    mv LAZYPIPE $OUTPUT
    ALIGN_TO_REF $OUTPUT/LAZYPIPE $OUTPUT/lazypipe-reconstructed.fa lazypipe
    cd $TOOL_PATH 
    conda activate base 

  fi
  #
  if [[ "$METASPADES" -eq "1" ]];
    then
    printf "Reconstructing with metaSPAdes\n\n"
    eval "$(conda shell.bash hook)"
    conda activate spades
    rm -rf metaspades_reconstruction
    mkdir metaspades_reconstruction
    cd metaspades_reconstruction

    mkdir input
    rm -rf output
    mkdir output	
    cp $READS1 $READS2 input
    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o metaspades-time.txt metaspades.py -o output -1 input/$aux_READS1 -2 input/$aux_READS2 -t $THREADS -m $MEMORY 
    mv metaspades-time.txt $OUTPUT
    mv output/scaffolds.fasta output/metaspades-reconstructed.fa
    cp output/metaspades-reconstructed.fa $OUTPUT
    cd ..
    conda activate base
    rm -rf metaspades_reconstruction
    cd $CURR_PATH
    CLASSIFY_INPUT metaSPAdes $OUTPUT/metaspades-reconstructed.fa
    mv metaSPAdes $OUTPUT
    ALIGN_TO_REF $OUTPUT/metaSPAdes $OUTPUT/metaspades-reconstructed.fa metaspades
    cd $TOOL_PATH 
  fi
  #
  if [[ "$METAVIRALSPADES" -eq "1" ]];
    then
    printf "Reconstructing with metaSPAdes\n\n"
    eval "$(conda shell.bash hook)"
    conda activate spades
    rm -rf metaviralspades_reconstruction
    mkdir metaviralspades_reconstruction
    cd metaviralspades_reconstruction

    mkdir input
    rm -rf output
    mkdir output	
    cp $READS1 $READS2 input
    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o metaviralspades-time.txt metaviralspades.py -o output -1 input/$aux_READS1 -2 input/$aux_READS2 -t $THREADS -m $MEMORY 
    mv metaviralspades-time.txt $OUTPUT
    mv output/scaffolds.fasta output/metaviralspades-reconstructed.fa
    cp output/metaviralspades-reconstructed.fa $OUTPUT
    cd ..
    conda activate base
    rm -rf metaviralspades_reconstruction
    cd $CURR_PATH
    CLASSIFY_INPUT metaviralSPAdes $OUTPUT/metaviralspades-reconstructed.fa
    mv metaviralSPAdes $OUTPUT
    ALIGN_TO_REF $OUTPUT/metaviralSPAdes $OUTPUT/metaviralspades-reconstructed.fa metaviralspades
    cd $TOOL_PATH
  fi
  #
  if [[ "$PEHAPLO" -eq "1" ]]; #err
    then
    printf "Reconstructing with PEHaplo\n\n"
    eval "$(conda shell.bash hook)"
    conda activate bio2
    cd TAR-VIR/PEHaplo/  
    rm -rf assembly
    mkdir assembly  
    cd assembly	    
    cp $READS1 $READS2 .
    sed -n '1~4s/^@/>/p;2~4p' $aux_READS1 > reads_1.fa
    sed -n '1~4s/^@/>/p;2~4p' $aux_READS2 > reads_2.fa      
    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o pehaplo-time.txt python ../pehaplo.py -f1 reads_1.fa -f2 reads_2.fa -l 10 -r 150 -t $THREADS -m ${MEMORY}GB -correct yes 
    cp pehaplo-time.txt $OUTPUT
    mv Contigs.fa pehaplo-reconstructed.fa
    cp pehaplo-reconstructed.fa $OUTPUT
    rm -rf Contigs.fa
    cd ../
    rm -rf assembly
    cd ../../
    conda activate base  
    cd $CURR_PATH
    CLASSIFY_INPUT PEhaplo $OUTPUT/pehaplo-reconstructed.fa
    mv PEhaplo $OUTPUT
    ALIGN_TO_REF $OUTPUT/PEhaplo $OUTPUT/pehaplo-reconstructed.fa pehaplo
    cd $TOOL_PATH
  fi
  #
  if [[ "$QURE" -eq "1" ]];
    then
    printf "Reconstructing with QuRe\n\n"
    eval "$(conda shell.bash hook)"
    conda activate java-env
    CREATE_PAIRED_FA_FILES
    cd QuRe_v0.99971/
    for virus in "${VIRUSES[@]}"
      do
      cp ../paired.fa $virus .
      aux_virus="$(cut -d'.' -f1 <<< $virus)"
      aux_virus="$(echo $aux_virus | awk -F/ '{print $NF}')"
      /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o qure-$aux_virus-time.txt java -Xmx${MEMORY}G -XX:MaxRAM=${MEMORY}G QuRe paired.fa $aux_virus.fa 
      mv paired_reconstructedVariants.txt results_$aux_virus.fa
    done  
    cat results_*.fa  >  qure-reconstructed.fa
    cp qure-reconstructed.fa $OUTPUT

    total_time=0
    total_mem=0
    total_cpu=0
    count=0

    for f in qure-*-time.txt
    do
      echo "Processing $f" 
      TIME=`cat $f | grep "TIME" | awk '{ print $2;}'`;
      MEM=`cat $f | grep "MEM" | awk '{ print $2;}'`;
      CPU=`cat $f | grep "CPU_perc" | awk '{ print $2;}'`;
      CPU="$(cut -d'%' -f1 <<< $CPU)"
      total_time=`echo "$total_time+$TIME" | bc -l`
      if [[ $MEM -gt $total_mem ]]
      then
        total_mem=$MEM
      fi
      total_cpu=`echo "$total_cpu+$CPU" | bc -l`
      count=`echo "$count+1" | bc -l`
    done
    printf "$total_cpu    -   $count     "
    total_cpu=$(echo $total_cpu \/ $count |bc -l | xargs printf %.0f)
    echo "TIME	$total_time
MEM	$total_mem
CPU_perc	$total_cpu%" > qure-time.txt
    mv qure-time.txt $OUTPUT
    rm -rf qure-*-time.txt
    rm paired_*
    rm *.fa
    cd ..
    conda activate base
    cd $CURR_PATH
    cp -r references QuRe
    mv QuRe $OUTPUT
    ALIGN_TO_REF $OUTPUT/QuRe $OUTPUT/qure-reconstructed.fa qure
    cd $TOOL_PATH
    fi
  #
  if [[ "$QVG" -eq "1" ]]; #unchanged
    then
    printf "Reconstructing with QVG\n\n"
    eval "$(conda shell.bash hook)"
    conda activate qvg
    cd QVG
    rm -rf reconstruction_files
    mkdir reconstruction_files
    rm -rf *.fa* 
    for virus in "${VIRUSES[@]}"
      do
      rm -rf dataset
      mkdir dataset
      echo "dataset" > dataset/samples
      cd dataset
      rm -rf output
      mkdir output
      cp $READS1 $READS2 .
      gzip -cvf $aux_READS1 > dataset_R1.fastq.gz
      gzip -cvf $aux_READS2 > dataset_R2.fastq.gz
      cd ..
      cp ${virus} reconstruction_files
      printf "$virus \n\n"
      aux_virus="$(cut -d'.' -f1 <<< $virus)"
      aux_virus="$(echo $aux_virus | awk -F/ '{print $NF}')"

      /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o qvg-$aux_virus-time.txt ./QVG.sh -r ./reconstruction_files/${aux_virus}.fa -samples-list ./dataset/samples -s ./dataset -o ./dataset/output -annot yes -np $THREADS
  
      cat dataset/output/samples_multifasta_* > dataset/output/qvg-${aux_virus}.fasta      
      cp dataset/output/qvg-${aux_virus}.fasta .
    done
    cat *.fasta > qvg-reconstructed.fa
    rm -rf *.fasta
    cp qvg-reconstructed.fa $OUTPUT 

    total_time=0
    total_mem=0
    total_cpu=0
    count=0
    for f in qvg-*-time.txt
    do
      echo "Processing $f" 
      TIME=`cat $f | grep "TIME" | awk '{ print $2;}'`;
      MEM=`cat $f | grep "MEM" | awk '{ print $2;}'`;
      CPU=`cat $f | grep "CPU_perc" | awk '{ print $2;}'`;
      CPU="$(cut -d'%' -f1 <<< $CPU)"
      total_time=`echo "$total_time+$TIME" | bc -l`
      if [[ $MEM -gt $total_mem ]]
      then
        total_mem=$MEM
      fi
      total_cpu=`echo "$total_cpu+$CPU" | bc -l`
      count=`echo "$count+1" | bc -l`

    printf "$total_cpu    -   $count     "
    total_cpu=$(echo $total_cpu \/ $count |bc -l | xargs printf %.0f)
    echo "TIME	$total_time
MEM	$total_mem
CPU_perc	$total_cpu%" > qvg-time.txt
    mv qvg-time.txt $OUTPUT
    rm -rf qvg-*-time.txt
    done
    conda activate base 
    cd ..
    cd $CURR_PATH
    cp -r references QVG
    mv QVG $OUTPUT
    ALIGN_TO_REF $OUTPUT/QVG $OUTPUT/qvg-reconstructed.fa qvg
    cd $TOOL_PATH
  fi
  #
  if [[ "$SPADES" -eq "1" ]];
    then
    printf "Reconstructing with SPAdes\n\n"
    eval "$(conda shell.bash hook)"
    conda activate spades
    rm -rf spades_reconstruction
    mkdir spades_reconstruction
    cd spades_reconstruction

    mkdir input
    rm -rf output
    mkdir output	
    cp $READS1 $READS2 input
    printf ".......   $aux_READS1 \n\n"
    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o spades-time.txt spades.py -o output -1 input/$aux_READS1 -2 input/$aux_READS2 -t $THREADS -m $MEMORY 
    mv spades-time.txt $OUTPUT
    mv output/scaffolds.fasta output/spades-reconstructed.fa
    cp output/spades-reconstructed.fa $OUTPUT
    cd ..
    conda activate base
    rm -rf spades_reconstruction
    #
    cd $CURR_PATH
    CLASSIFY_INPUT SPAdes $OUTPUT/spades-reconstructed.fa
    mv SPAdes $OUTPUT
    ALIGN_TO_REF $OUTPUT/SPAdes $OUTPUT/spades-reconstructed.fa spades
    cd $TOOL_PATH   
  fi
  #
  if [[ "$SSAKE" -eq "1" ]];
    then
    printf "Reconstructing with SSAKE\n\n"
    cd ssake/tools/
    cp $READS1 $READS2 .
    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o ssake-time.txt ./runSSAKE.sh $aux_READS1 $aux_READS2 10 assembly
    mv assembly_scaffolds.fa ssake-reconstructed.fa
    cp ssake-reconstructed.fa $OUTPUT
    mv ssake-time.txt $OUTPUT
    rm *.fa
    rm *.fq
    rm assembly*
    cd ../../
    cd $CURR_PATH
    CLASSIFY_INPUT SSAKE $OUTPUT/ssake-reconstructed.fa
    mv SSAKE $OUTPUT
    ALIGN_TO_REF $OUTPUT/SSAKE $OUTPUT/ssake-reconstructed.fa ssake
    cd $TOOL_PATH 
  fi
  #
  if [[ "$TRACESPIPE" -eq "1" ]];
    then
    printf "Reconstructing with TRACESPipe\n\n"
    eval "$(conda shell.bash hook)"
    conda activate tracespipe
    lzma -d VDB.fa.lzma
    cd tracespipe/
    cd input_data/
    cp $READS1 $READS2 .
    rm -rf *.fq.gz
    gzip $aux_READS1
    gzip $aux_READS2
    cd ../meta_data/
    echo "x:$aux_READS1.gz:$aux_READS2.gz" > meta_info.txt
    cd ../src/
    cp ../../VDB.fa .

    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o tracespipe-time.txt ./TRACESPipe.sh --run-meta --run-all-v-alig --very-sensitive -t $THREADS --cache 10
    cp tracespipe-time.txt ../
    cd ../output_data/TRACES_viral_consensus
    cat *.fa > ../../tracespipe-reconstructed.fa     
    cd ../../

    mv tracespipe-reconstructed.fa $OUTPUT
    mv tracespipe-time.txt $OUTPUT
    rm *-time.txt
    cd output_data
    mv TRACES_results/REPORT_META_VIRAL_ALL_SAMPLES.txt $OUTPUT    
    rm -rf *
    cd ..
    cd meta_data
    rm meta_info.txt
    cd ../input_data
    rm *
    cd ../../   
    CLASSIFY_TRACES $OUTPUT/REPORT_META_VIRAL_ALL_SAMPLES.txt $OUTPUT/TRACESPipe
    ALIGN_TO_REF $OUTPUT/TRACESPipe $OUTPUT/tracespipe-reconstructed.fa tracespipe
    cd $TOOL_PATH
    conda activate base 
  fi
  #
  if [[ "$TRACESPIPELITE" -eq "1" ]];
    then
    printf "Reconstructing with TRACESPipeLite\n\n"
    eval "$(conda shell.bash hook)"  
    conda activate tracespipelite
    cd TRACESPipeLite/src/  	
    cp $READS1 $READS2 .
    lzma -d VDB.mfa.lzma
    rm -rf *.fq.gz
    gzip $aux_READS1
    gzip $aux_READS2
    /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o tracespipelite-time.txt ./TRACESPipeLite.sh --similarity 5 --threads $THREADS --reads1 $aux_READS1.gz --reads2 $aux_READS2.gz --database VDB.mfa --output test_viral_analysis --no-plots
    cd test_viral_analysis
    for virus in $(ls)
    do
      if [ -d $virus ] 
      then
        printf "copying $virus\n\n"
        cd $virus*   
        cp *-consensus.fa ../../
        cd ..
      fi 
    done
    cd ..
    cat *-consensus.fa > tracespipelite-reconstructed.fa
    mv tracespipelite-time.txt $OUTPUT
    cp tracespipelite-reconstructed.fa $OUTPUT
    rm *-*.fa
    mv test_viral_analysis/best-viral-metagenomics.txt $OUTPUT
    rm -rf test_viral_analysis*
    rm *.fq.gz
    rm *-consensus.fa
    cd ../../
    conda activate base
    CLASSIFY_TRACES $OUTPUT/best-viral-metagenomics.txt $OUTPUT/TRACESPipeLite
    ALIGN_TO_REF $OUTPUT/TRACESPipeLite $OUTPUT/tracespipelite-reconstructed.fa tracespipelite
    cd $TOOL_PATH
  fi
  #
  if [[ "$VIRGENA" -eq "1" ]]; #unchanged
    then
  printf "Reconstructing with VirGenA\n\n"
  eval "$(conda shell.bash hook)"
  conda activate java-env
  cd release_v1.4
  chmod +x ./tools/vsearch
    
    for virus in "${VIRUSES[@]}"
    do
      #rm -rf ${dataset}_*.fq
      rm -rf ${virus}.fa
      rm -rf *.gz
      
      cp $READS1 $READS2 .
      cp ${virus} .
    
      gzip *.fq
      
      #rm -rf $dataset-$virus 
      #mkdir $dataset-$virus
      
      aux_virus="$(cut -d'.' -f1 <<< $virus)"
      aux_virus="$(echo $aux_virus | awk -F/ '{print $NF}')"
    
      echo "<config>
    <Data>
        <pathToReads1>$aux_READS1.gz</pathToReads1>
        <pathToReads2>$aux_READS2.gz</pathToReads2>
        <InsertionLength>1000</InsertionLength>
    </Data>
    <Reference>${aux_virus}.fa</Reference>
    <OutPath>./res/$aux_virus</OutPath>
    <ThreadNumber>-1</ThreadNumber>
	<BatchSize>1000</BatchSize>
    <ReferenceSelector>
		<Enabled>true</Enabled>
        <UseMajor>false</UseMajor>
        <ReferenceMSA>./data/HIV_RefSet_msa.fasta</ReferenceMSA>
        <PathToVsearch>./tools/vsearch</PathToVsearch>
        <UclustIdentity>0.95</UclustIdentity>
        <MinReadLength>50</MinReadLength>
        <MinContigLength>1000</MinContigLength>
        <Delta>0.05</Delta>
		<MaxNongreedyComponentNumber>5</MaxNongreedyComponentNumber>
        <MapperToMSA>
			<K>7</K>
			<IndelToleranceThreshold>1.5</IndelToleranceThreshold>
			<pValue>0.01</pValue>
			<RandomModelParameters>
				<Order>4</Order>
				<ReadNum>1000</ReadNum>
				<Step>10</Step>
			</RandomModelParameters>
        </MapperToMSA>
        <Graph>
            <MinReadNumber>5</MinReadNumber>
            <VertexWeight>10</VertexWeight>
			<SimilarityThreshold>0.05</SimilarityThreshold>
			<Debug>false</Debug>
        </Graph>
        <Debug>false</Debug>
    </ReferenceSelector>
    <Mapper>
		<K>5</K>
		<IndelToleranceThreshold>1.25</IndelToleranceThreshold>
        <pValue>0.01</pValue>
		<RandomModelParameters>
			<Order>4</Order>
			<ReadNum>1000</ReadNum>
			<Step>10</Step>
		</RandomModelParameters>
        <Aligner>
            <Match>2</Match>
            <Mismatch>-3</Mismatch>
            <GapOpenPenalty>5</GapOpenPenalty>
            <GapExtensionPenalty>2</GapExtensionPenalty>
        </Aligner>
    </Mapper>
    <ConsensusBuilder>
        <IdentityThreshold>0.9</IdentityThreshold>
        <CoverageThreshold>0</CoverageThreshold>
        <MinIntersectionLength>10</MinIntersectionLength>
        <MinTerminationReadsNumber>1</MinTerminationReadsNumber>
        <Reassembler>
            <IdentityThreshold>0.9</IdentityThreshold>
            <MinTerminatingSequenceCoverage>0</MinTerminatingSequenceCoverage>
            <PairReadTerminationThreshold>0.1</PairReadTerminationThreshold>
            <MinReadLength>50</MinReadLength>
        </Reassembler>
        <Debug>false</Debug>
    </ConsensusBuilder>
    <Postprocessor>
        <Enabled>true</Enabled>
        <MinFragmentLength>500</MinFragmentLength>
        <MinIdentity>0.99</MinIdentity>
		<MinFragmentCoverage>0.99</MinFragmentCoverage>
        <Debug>false</Debug>
    </Postprocessor>
</config>" > conf.xml
    
    
    
    timeout --signal=SIGINT ${VIRGENA_TIMEOUT}m /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o virgena-$aux_virus-time.txt java -jar VirGenA.jar assemble -c conf.xml # config_test_linux.xml
    #java -jar VirGenA.jar map -c config.xml -r ../B19.fa -p1 ../DS1_1.fq -p2 ../DS1_2.fq
    #rm $virus*
    
    done
    cd res
    cat *_complete_genome_assembly.fasta > virgena-reconstructed.fa
    rm -rf *_complete_genome_assembly.fasta
    
    
    total_time=0
    total_mem=0
    total_cpu=0 
    count=0
    for f in ../virgena-*-time.txt
    do
      echo "Processing $f" 
      TIME=`cat $f | grep "TIME" | awk '{ print $2;}'`;
      MEM=`cat $f | grep "MEM" | awk '{ print $2;}'`;
      CPU=`cat $f | grep "CPU_perc" | awk '{ print $2;}'`;
      CPU="$(cut -d'%' -f1 <<< $CPU)"       
      total_time=`echo "$total_time+$TIME" | bc -l`      
      if [[ $MEM -gt $total_mem ]]
      then
        total_mem=$MEM
      fi
      total_cpu=`echo "$total_cpu+$CPU" | bc -l`
      count=`echo "$count+1" | bc -l`     
    done
    printf "$total_cpu    -   $count     "
    total_cpu=$(echo $total_cpu \/ $count |bc -l | xargs printf %.0f)
    echo "TIME	$total_time
MEM	$total_mem
CPU_perc	$total_cpu%" > ../virgena-time.txt
        
    mv ../virgena-time.txt $OUTPUT
    mv virgena-reconstructed.fa $OUTPUT
    rm -rf *
    cd ..
    rm virgena-time.txt
    
  cd ..
  conda activate base
  
    cd $CURR_PATH
    CLASSIFY_INPUT VirGenA $OUTPUT/virgena-reconstructed.fa
    mv VirGenA $OUTPUT
    ALIGN_TO_REF $OUTPUT/VirGenA $OUTPUT/virgena-reconstructed.fa virgena
    cd $TOOL_PATH
  fi
  #
  if [[ "$VISPA" -eq "1" ]];
    then
    printf "Reconstructing with ViSpA\n\n"  
    eval "$(conda shell.bash hook)"
    CREATE_PAIRED_FA_FILES
    conda activate vispa  
    cd home
    rm -rf test
    mkdir test
    cd test 	
    cp ../../paired.fa .
    for virus in "${VIRUSES[@]}"
    do
      cp $virus .
      #read a
      aux_virus="$(cut -d'.' -f1 <<< $virus)"
      aux_virus="$(echo $aux_virus | awk -F/ '{print $NF}')"

      echo "" >> data.txt

      /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o vispa-$aux_virus-time.txt  ../code/vispa_mosaik/./main_mosaik.bash paired.fa $aux_virus.fa $THREADS 100 120
      mv paired_I_*_*_CNTGS_DIST0.txt tmp_$aux_virus.fa
      done

      for f in tmp_*.fa
        do
        printf "changing $f \n\n"
        content=$(cat $f)
        echo ">${f}
${content}" > zz_$f
        done

      cat zz_tmp_*.fa > vispa-reconstructed.fa
      cp vispa-reconstructed.fa $OUTPUT

      total_time=0
      total_mem=0
      total_cpu=0 
      count=0
      for f in vispa-*-time.txt
      do
        echo "Processing $f" 
        TIME=`cat $f | grep "TIME" | awk '{ print $2;}'`;
        MEM=`cat $f | grep "MEM" | awk '{ print $2;}'`;
        CPU=`cat $f | grep "CPU_perc" | awk '{ print $2;}'`;
        CPU="$(cut -d'%' -f1 <<< $CPU)"       
        total_time=`echo "$total_time+$TIME" | bc -l`      
        if [[ $MEM -gt $total_mem ]]
        then
          total_mem=`echo "$total_mem+$MEM" | bc -l`
        fi
        total_cpu=`echo "$total_cpu+$CPU" | bc -l`
        count=`echo "$count+1" | bc -l`     
      done
    total_cpu=$(echo $total_cpu \/ $count |bc -l | xargs printf %.0f)
    echo "TIME	$total_time
MEM	$total_mem
CPU_perc	$total_cpu%" > vispa-time.txt
    cp vispa-time.txt $OUTPUT
    rm -rf vispa-*-time.txt
    rm *.fa
    cd ../
    rm *.fa
    rm *.fq.gz
    cd res
    rm *.fa*
    rm new*.txt
    cd ../../

    conda activate base  
    cd $CURR_PATH
    cp -r references ViSpA
    mv ViSpA $OUTPUT
    ALIGN_TO_REF $OUTPUT/ViSpA $OUTPUT/vispa-reconstructed.fa vispa
    cd $TOOL_PATH
  fi
  #
  if [[ "$VPIPE" -eq "1" ]];
    then
    printf "Reconstructing with V-pipe\n\n"
    eval "$(conda shell.bash hook)"
    conda activate V-pipe
    cd V-pipe/
    mkdir references 
    for virus in "${VIRUSES[@]}"
      do 
      aux_virus="$(cut -d'.' -f1 <<< $virus)"
      aux_virus="$(echo $aux_virus | awk -F/ '{print $NF}')"
      cd config/
      echo "---
name: $aux_virus
general:
    aligner: "bwa"
input:
    reference: \"{VPIPE_BASEDIR}/../resources/$aux_virus/$aux_virus.fa\"
snv:
    consensus: false
lofreq:
    consensus: true" > $aux_virus.yaml
      cd ../resources/
      rm -rf $aux_virus
      mkdir $aux_virus
      cp $virus $aux_virus
      cd ../references/
      cp $virus .
      cd ..
      echo "general:
  virus_base_config: '$aux_virus'
  
input:
  datadir: samples
  samples_file: $(pwd)/config/samples.tsv
  reference: $(pwd)/references/${aux_virus}.fa
  trim_percent_cutoff: 0.1
output:
  datadir: $(pwd)/results
  snv: false
  local: false
  global: false
  visualization: false
  QA: false" > config/config.yaml 
      echo "SRR10903401	20200102" > config/samples.tsv
      rm -rf samples
      mkdir samples
      cd samples
      mkdir SRR10903401
      cd SRR10903401
      mkdir 20200102
      cd 20200102
      mkdir raw_data
      cd raw_data 
      cp $READS1 $READS2 .
      sed -i 's/J/I/g' $aux_READS1  
      sed -i 's/J/I/g' $aux_READS2  

      mv $aux_READS1 SRR10903401_R1.fq
      mv $aux_READS2 SRR10903401_R2.fq
      cd ../../../

      cd ..

      /bin/time -f "TIME\t%e\nMEM\t%M\nCPU_perc\t%P" -o v-pipe-$aux_virus-time.txt ./vpipe  --cores $THREADS --conda-frontend conda
      cp results/SRR10903401/20200102/references/ref_majority.fasta .
      mv ref_majority.fasta $aux_virus.fasta

      cd references
      rm ${aux_virus}*
      cd ..
      cd resources
      rm -rf ${aux_virus}
      cd results
      rm -rf *
      cd ..
      cd samples
      rm -rf *
      cd ..
      cd config
      rm ${aux_virus}.yaml
      cd ..
      
      
    done

    total_time=0
    total_mem=0
    total_cpu=0
    count=0
    for f in v-pipe-*-time.txt
    do
      echo "Processing $f" 
      TIME=`cat $f | grep "TIME" | awk '{ print $2;}'`;
      MEM=`cat $f | grep "MEM" | awk '{ print $2;}'`;
      CPU=`cat $f | grep "CPU_perc" | awk '{ print $2;}'`;
      CPU="$(cut -d'%' -f1 <<< $CPU)"
      total_time=`echo "$total_time+$TIME" | bc -l`
      if [[ $MEM -gt $total_mem ]]
      then
        total_mem=$MEM
      fi
      total_cpu=`echo "$total_cpu+$CPU" | bc -l`
      count=`echo "$count+1" | bc -l`
    done
    printf "$total_cpu    -   $count     "
    total_cpu=$(echo $total_cpu \/ $count |bc -l | xargs printf %.0f)
    echo "TIME	$total_time
MEM	$total_mem
CPU_perc	$total_cpu%" > v-pipe-time.txt

    mv v-pipe-time.txt $OUTPUT 
    rm -rf v-pipe-*-time.txt
    cat *.fasta > v-pipe-reconstructed.fa
    rm -rf *.fasta
    
    cat v-pipe-reconstructed.fa | tr [:lower:] [:upper:] > tmp.txt
    mv tmp.txt v-pipe-reconstructed.fa
      
    mv v-pipe-reconstructed.fa $OUTPUT    

    cd ../
    conda activate base
    cd $CURR_PATH
    cp -r references V-pipe
    mv V-pipe $OUTPUT
    ALIGN_TO_REF $OUTPUT/V-pipe $OUTPUT/v-pipe-reconstructed.fa v-pipe
    cd $TOOL_PATH
  fi 
  cd $CURR_PATH
  PERFORM_MULTIPLE_ALIGNMENT $OUTPUT/consensus
  CREATE_FINAL_CONSENSUS $OUTPUT/consensus
fi
#
################################################################################
#
#PERFORM_MULTIPLE_ALIGNMENT $OUTPUT/consensus
if [[ "$AWK" -eq "1" ]];
  then
  CREATE_FINAL_CONSENSUS $OUTPUT/consensus
fi
#
if [[ "$EMBOSS" -eq "1" ]];
  then
  CREATE_FINAL_CONSENSUS $OUTPUT/consensus
fi
#
#printf "Evaluating results based on the classification made of the input reads.\n\n"
#cd references
#cat *.fa > reference_file.fa
#cp reference_file.fa ../ 
#cd ..
#./Evaluation.sh -ref reference_file.fa
#
################################################################################
#

