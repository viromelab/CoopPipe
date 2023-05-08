#!/bin/bash
#
#
D_PATH="out_analysis"
REF_FILE=""
#
declare -a ANALYSIS=("tracespipelite" "spades" "metaspades" "metaviralspades" "coronaspades" "ssake" "tracespipe" "lazypipe" "pehaplo" "haploflow");
declare -a NO_ANALYSIS=("qvg" "qure" "vispa" "virgena" "v");
#
declare -a CLASSIFICATION=("tracespipelite" "tracespipe" "lazypipe");
declare -a NO_CLASSIFICATION=("spades" "metaspades" "metaviralspades" "coronaspades" "ssake" "pehaplo" "qvg" "qure" "vispa" "virgena" "haploflow" "v");
#
################################################################################
#
SHOW_MENU () {
  echo " --------------------------------------------------------- ";
  echo "                                                           ";
  echo " Evaluation.sh : TRACESPipe lite version v2.1              ";
  echo "                                                           ";
  echo " This is a lite version of TRACESPipe. It provides         ";
  echo " automatic reconstruction (reference-based only) of        ";
  echo " viral genomes and performs basic analyses.                ";
  echo "                                                           ";
  echo " Program options ----------------------------------------- ";
  echo "                                                           ";
  echo " -h, --help                   Show this,                   ";
  echo " -ref <STR>, --reference <STR>         Reference for the   ";
  echo "                                      reconstruction       ";
  echo " -dir<STR>, --dir_reconstructed <STR>  Directory where the ";
  echo "                  results of the reconstruction are stored,";
  echo " Example ------------------------------------------------  ";
  echo "                                                           ";
  echo " ./Evaluation.sh -ref reference.fa -dir out_analysis       ";
  echo "                                                           ";
  echo " --------------------------------------------------------  ";
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
    -ref|--reference)
      REF_FILE="$2";
      shift 2;
    ;;
    -dir|--dir_reconstructed)
      D_PATH="$2";
      shift 2;
    ;;
    -*) # unknown option with small
    echo "Invalid arg ($1)!";
    echo "For help, try: TRACESPipeLite.sh -h"
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
if [ -f "$REF_FILE" ];
  then
  echo "File	Time(s)	SNPs	AvgIdentity	NCD	NRC	Mem(GB)	%CPU	Nr contigs	Metagenomic_analysis	Metagenomic_classification	Coverage	SNP_mutations_DS	Contamination_ds" > total_stats.tsv
#
  REF_FILE="$(pwd)/$(echo $REF_FILE | awk -F/ '{print $NF}')"
  
  cd $D_PATH
  for file in `ls -1 *-reconstructed.fa*` #for each fasta file in curr dir
  do 
    rm -rf out.report	 
    TIME=-1
    SNPS=-1
    IDEN=1
    NCD=1
    NRC=1
    MEM=-1
    CPU_P=-1
    NR_SPECIES=-1
    DOES_ANALYSIS=-1
    DOES_CLASSIFICATION=-1
  
    fst_char=$(cat $file | head -c 1)
    if [[ -z "$fst_char" ]]; then
      printf "The result file is empty."    
    else
      dos2unix $file  
      awk -i inplace '{ while(sub(/QuRe./,int(rand()*99999999999)+1)); print }' $file
      awk -i inplace '{ while(sub(/results/,int(rand()*99999999999)+1)); print }' $file

      dnadiff $file $REF_FILE; #run dnadiff
      
      IDEN=`cat out.report | grep "AvgIdentity " | head -n 1 | awk '{ print $2;}'`;  #retrieve results
      ALBA=`cat out.report | grep "AlignedBases " | head -n 1 | awk '{ print $2;}'`;
      SNPS=`cat out.report | grep TotalSNPs | awk '{ print $2;}'`;
      TBASES=`cat out.report | grep "TotalBases" | awk '{ print $2;}'`;
      AUX="$(cut -d')' -f1 <<< "$ALBA")"
      PERALBA="$(cut -d'(' -f2 <<< "$AUX")"
      TALBA="$(cut -d'(' -f1 <<< "$ALBA")" 
         
      NRBASES=`cat out.report | grep "TotalBases" | awk '{ print $2;}'`;  
          
      file_wout_extension="$(cut -d'.' -f1 <<< $file)"
      file_wout_extension="$(cut -d'-' -f1 <<< $file_wout_extension)"
      
      if [ "$file_wout_extension" == "v" ];
        then
        file_wout_extension="v-pipe"
      fi
      
      TIME=`cat $file_wout_extension-time.txt | grep "TIME" | awk '{ print $2;}'`;
      MEM=`cat $file_wout_extension-time.txt | grep "MEM" | awk '{ print $2;}'`;
      CPU_P=`cat $file_wout_extension-time.txt | grep "CPU_perc" | awk '{ print $2;}'`;
      #TMP=$(($TALBA * 100))
      #ACCURACY=$(echo $TMP \/ $TBASES |bc -l | xargs printf %.3f)
      
      NAME_TOOL="$(cut -d'-' -f1 <<< $file_wout_extension)"
      
      DOES_ANALYSIS="?"
      DOES_CLASSIFICATION="?"      

      for i in "${ANALYSIS[@]}" #check if the tool does metagenomic analysis
      do
        if [ "$i" == "$NAME_TOOL" ] ; then
          DOES_ANALYSIS="Yes"
          break
        fi 
      done
      
      if [ "$DOES_ANALYSIS" == "?" ] ; then
        for i in "${NO_ANALYSIS[@]}"
        do
          if [ "$i" == "$NAME_TOOL" ] ; then
            DOES_ANALYSIS="No"
            break
          fi 
        done
      fi
      
      for i in "${CLASSIFICATION[@]}" #check if the tool does metagenomic classification
      do
        if [ "$i" == "$NAME_TOOL" ] ; then
          DOES_CLASSIFICATION="Yes"
          break
        fi
      done
      
      if [ "$DOES_CLASSIFICATION" == "?" ] ; then
        for i in "${NO_CLASSIFICATION[@]}"
        do
          if [ "$i" == "$NAME_TOOL" ] ; then
            DOES_CLASSIFICATION="No"
            break
          fi 
        done
      fi
      
      NR_SPECIES=$(grep '>' $file -c)
      gto_fasta_rand_extra_chars < $file > tmp.fa
      gto_fasta_to_seq < tmp.fa > $file_wout_extension-reconstructed.seq
      
      #Compressing sequences C(X) or C(X,Y)
      GeCo3 -tm 1:1:0:1:0.9/0:0:0 -tm 7:10:0:1:0/0:0:0 -tm 16:100:1:10:0/3:10:0.9 -lr 0.03 -hs 64 $file_wout_extension-reconstructed.seq 
      COMPRESSED_SIZE_WOUT_REF=$(ls -l $file_wout_extension-reconstructed.seq.co | cut -d' ' -f5)
      rm $file_wout_extension*.seq.*
      
      #Conditional compression C(X|Y) [use reference and target]

      GeCo3 -rm 20:500:1:12:0.9/3:100:0.9 -rm 13:200:1:1:0.9/0:0:0 -tm 1:1:0:1:0.9/0:0:0 -tm 7:10:0:1:0/0:0:0 -tm 16:100:1:10:0/3:10:0.9 -lr 0.03 -hs 64 -r $REF_FILE $file_wout_extension-reconstructed.seq
      COMPRESSED_SIZE_COND_COMPRESSION=$(ls -l $file_wout_extension-reconstructed.seq.co | cut -d' ' -f5)   
      rm $file_wout_extension*.seq.*
      
      #Relative compression (only reference models) C(X||Y)
      GeCo3 -rm 20:500:1:12:0.9/3:100:0.9 -rm 13:200:1:1:0.9/0:0:0 -lr 0.03 -hs 64 -r $REF_FILE $file_wout_extension-reconstructed.seq
      COMPRESSED_SIZE_W_REF=$(ls -l $file_wout_extension-reconstructed.seq.co | cut -d' ' -f5)      
      rm $file_wout_extension*.seq.*            
      FILE_SIZE=$(ls -l $file_wout_extension-reconstructed.fa | cut -d' ' -f5)
     
      NCD=$(echo $COMPRESSED_SIZE_COND_COMPRESSION \/ $COMPRESSED_SIZE_WOUT_REF |bc -l | xargs printf %.3f)
       
      AUX_MULT=$(echo "$FILE_SIZE * 2" | bc -l )
      NRC=$(echo $COMPRESSED_SIZE_W_REF \/ $AUX_MULT|bc -l | xargs printf %.3f)      
      
      IDEN=$(echo $IDEN |bc -l | xargs printf %.3f)
      MEM=$(echo $MEM \/ 1048576 |bc -l | xargs printf %.3f)
      
    #file	exec_time	snps	avg_identity	NCD	NRC	max_mem	cpu_avg	nr_contigs_reconstructed	metagenomic_analysis	metagenomic_classification	coverage	snp_dataset
    echo "$file	$TIME	$SNPS	$IDEN	$NCD	$NRC	$MEM	$CPU_P	$NR_SPECIES	$DOES_ANALYSIS	$DOES_CLASSIFICATION	$coverage	$snp_ds	$cnt_ds" >> ../total_stats.tsv   
    fi 
done
cd ..
else
  printf "Reference is invalid!\n\n"
fi
#
