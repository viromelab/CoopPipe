#!/bin/bash
#
declare -a VIRUSES_AVAILABLE=("B19V" "BuV" "CuV" "HBoV" "AAV" "BKPyV" "JCPyV" "KIPyV"
                    "WUPyV" "MCPyV" "HPyV6" "HPyV7" "TSPyV" "HPyV9" "MWPyV"
                    "STLPyV" "HPyV12" "NJPyV" "LIPyV" "SV40" "TTV" "TTVmid"
                    "TTVmin" "HAV" "HBV" "HCV" "HDV" "HEV" "SENV" "HPV11" "HPV16" 
                    "HPV18" "HPV31" "HPV39" "HPV45" "HPV51" "HPV56" "HPV58" "HPV59"
                    "HPV68" "HPV77""HPV2" "HPV6"  "HSV-1"
                    "HSV-2" "VZV" "EBV" "HCMV" "HHV6" "HHV7" "KSHV" "ReDoV"
                    "VARV" "MPXV" "EV" "SARS2" "HERV" "MT");
#
rm -rf Graphs_k
mkdir Graphs_k
#
nr_virus=0;
#
cd Graphs_k
cp ../Results/total_stats_k.tsv .
#
for virus in "${VIRUSES_AVAILABLE[@]}" 
  do
  
  sort -k3 -g total_stats_k.tsv > tmp.tsv
  data=$(cat tmp.tsv | tr ',' '.' | grep -w "${virus}" > "$virus")
  
  if [ -s "${virus}" ];
    then 
    printf "$virus\n\n"
    
    printf "1\n"
    gnuplot << EOF
      reset
      set terminal pdfcairo enhanced color font 'Verdade,9'
      set output "Identity_$virus.pdf"
      set datafile separator "\t"
    
      set autoscale x
      set autoscale y
      
      set offsets graph 0, 0, 0.05, 0.05
   
      set ytics auto
      set xtics auto
      set ylabel "Average Identity"
      set xlabel "K"
      set multiplot layout 1,1
      set rmargin 5
      set key at screen 1, graph 1 
      
      plot "$virus" u 3:5 with linespoints lc 0 notitle
  
EOF
    printf "2\n"
    gnuplot << EOF
      reset
      set terminal pdfcairo enhanced color font 'Verdade,9'
      set output "NCSD_$virus.pdf"
      set datafile separator "\t"
    
      set autoscale x
      set autoscale y
      
      set offsets graph 0, 0, 0.05, 0.05
   
      set ytics auto
      set xtics auto
      set ylabel "NCSD"
      set xlabel "K"
      set multiplot layout 1,1
      set rmargin 5
      set key at screen 1, graph 1 
        
      plot "$virus" u 3:6 with linespoints lc 0 notitle
  
EOF
    printf "3\n"
    gnuplot << EOF
      reset
      set terminal pdfcairo enhanced color font 'Verdade,9'
      set output "NRC_$virus.pdf"
      set datafile separator "\t"
    
      set autoscale x
      set autoscale y
      
      set offsets graph 0, 0, 0.05, 0.05
   
      set ytics auto
      set xtics auto
      set ylabel "NRC"
      set xlabel "K"
      set multiplot layout 1,1
      set rmargin 5
      set key at screen 1, graph 1 
        
      plot "$virus" u 3:7 with linespoints lc 0 notitle
  
EOF
  
  fi
  
  rm $virus
#
done
#
#rm *.tsv
cd ..
#
