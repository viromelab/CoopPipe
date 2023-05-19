#!/bin/bash
#
declare -a VIRUSES_AVAILABLE=("B19V" "BuV" "CuV" "HBoV" "AAV" "BKPyV" "JCPyV" "KIPyV"
                    "WUPyV" "MCPyV" "HPyV6" "HPyV7" "TSPyV" "HPyV9" "MWPyV"
                    "STLPyV" "HPyV12" "NJPyV" "LIPyV" "SV40" "TTV" "TTVmid"
                    "TTVmin" "HAV" "HBV" "HCV" "HDV" "HEV" "SENV" "HPV2"
                    "HPV6" "HPV11" "HPV16" "HPV18" "HPV31" "HPV39" "HPV45"
                    "HPV51" "HPV56" "HPV58" "HPV59" "HPV68" "HPV77" "HSV-1"
                    "HSV-2" "VZV" "EBV" "HCMV" "HHV6" "HHV7" "KSHV" "ReDoV"
                    "VARV" "MPXV" "EV" "SARS2" "HERV" "MT");
#
rm -rf Graphs
mkdir Graphs
#
rm -rf avg_identity
rm -rf NCD
rm -rf NRC
mkdir avg_identity
mkdir NCD
mkdir NRC
#
nr_virus=0;
#
for virus in "${VIRUSES_AVAILABLE[@]}" 
  do
  
  
  cd avg_identity
  cp ../Results/total_stats.tsv .
 
  sort -k5 -r -n total_stats.tsv > tmp.tsv
  cooppipe=$(cat tmp.tsv | tr ',' '.' | grep "${virus}" | grep -w "cooppipealt-$virus-consensus.fa" > "${virus}")
  best_tool=$(cat tmp.tsv | tr ',' '.'  | grep "${virus}" | grep -w --invert-match "cooppipealt-$virus-consensus.fa" | head -1 >> "${virus}")
  
  if [ -s ${virus} ];
    then
    nr_virus=$(echo $nr_virus + 1 |bc -l)
  fi
  
  cd ..
  
  cd NRC
  cp ../Results/total_stats.tsv .
 
  sort -k6 -n total_stats.tsv > tmp.tsv
  cooppipe=$(cat tmp.tsv | tr ',' '.' | grep "${virus}" | grep -w "cooppipealt-$virus-consensus.fa" > "${virus}")
  best_tool=$(cat tmp.tsv | tail -n +2 | tr ',' '.'  | grep "${virus}" | grep -w --invert-match "cooppipealt-$virus-consensus.fa" | head -1 >> "${virus}")
  
  cd ..
  
  cd NCD
  cp ../Results/total_stats.tsv .
 
  sort -k7 -n total_stats.tsv > tmp.tsv
  cooppipe=$(cat tmp.tsv | tr ',' '.' | grep "${virus}" | grep -w "cooppipealt-$virus-consensus.fa" > "${virus}")
  best_tool=$(cat tmp.tsv | tail -n +2 | tr ',' '.'  | grep "${virus}" | grep -w --invert-match "cooppipealt-$virus-consensus.fa" | head -1 >> "${virus}")
 
  cd ..  
done
#
cd avg_identity
find . -type f -empty -delete
rm *.tsv
cd ..
#
cd NCD
find . -type f -empty -delete
rm *.tsv
cd ..
#
cd NRC
find . -type f -empty -delete
rm *.tsv
cd ..
#
list_avg=($(ls avg_identity))  
list_ncd=($(ls NCD))
list_nrc=($(ls NRC))
#
#  
#printf "${list[*]} \n\n"
#
cd avg_identity
#
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Verdade,9'
    set output "Identity.pdf"
    set datafile separator "\t"
    
    spacing_x = 30
    ymax = 100.1
    ymin = 98.5
    offset = ( ymax - ymin ) / 15.0    
    set yrange [ymin:ymax]
    set xrange [0:$nr_virus * spacing_x]
   
    set ytics auto
    set ylabel "Average Identity"
    set xlabel "Viruses"
    set multiplot layout 1,1
    set rmargin 5
    set key at screen 1, graph 1  
    
    count = 10
    do for [ file in "${list_avg[@]}"]{  
      set xtics (file count)
      plot file using (count):5:13 with labels point  pt 7 offset char 6,-0.1 notitle
      count = count + spacing_x
      ymax = ymax - offset
    }
EOF
#
cp *.pdf ../Graphs
cd ../NCD
#
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Verdade,9'
    set output "NCD.pdf"
    set datafile separator "\t"
    
    spacing_x = 30
    ymax = 0.3
    ymin = 0
    offset = ( ymax - ymin ) / 15.0    
    set yrange [ymin:ymax]
    set xrange [0:$nr_virus * spacing_x]
   
    set ytics auto
    set ylabel "NCD"
    set xlabel "Viruses"
    set multiplot layout 1,1
    set rmargin 5
    set key at screen 1, graph 1  
    
    count = 10
    do for [ file in "${list_ncd[@]}"]{  
      set xtics (file count)
      plot file using (count):6:13 with labels point  pt 7 offset char 6,-0.1 notitle
      count = count + spacing_x
    }
EOF
#
cp *.pdf ../Graphs
cd ../NRC
#
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Verdade,9'
    set output "NRC.pdf"
    set datafile separator "\t"
    
    spacing_x = 30
    ymax = 0.15
    ymin = 0
    offset = ( ymax - ymin ) / 15.0    
    set yrange [ymin:ymax]
    set xrange [0:$nr_virus * spacing_x]
   
    set ytics auto
    set ylabel "NCD"
    set xlabel "Viruses"
    set multiplot layout 1,1
    set rmargin 5
    set key at screen 1, graph 1  
    
    count = 10
    do for [ file in "${list_nrc[@]}"]{  
      set xtics (file count)
      plot file using (count):7:13 with labels point  pt 7 offset char 6,-0.1 notitle
      count = count + spacing_x
    }
EOF
#
cp *.pdf ../Graphs
cd ..
#
printf "Finished generating graphs!\n"
#
#set multiplot
#plot file u 13:5 title file with linespoints linestyle count
#rm -rf total_stats.tsv
#, "" u 13:5:5 with labels offset char 0,1
#plot "total_stats.tsv" using 13:6:(sprintf("(%d, %d)", $1, $2)) with labels notitle
#set key at 100., 100.
