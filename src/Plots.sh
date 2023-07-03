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
rm -rf Graphs
mkdir Graphs
#
rm -rf avg_identity
rm -rf NCSD
rm -rf NRC
mkdir avg_identity
mkdir NCSD
mkdir NRC
#
nr_virus=0;
#
for virus in "${VIRUSES_AVAILABLE[@]}" 
  do
  
  
  cd avg_identity
  cp ../Results/total_stats.tsv .
 
  sort -k5 -r -g total_stats.tsv > tmp.tsv
  cooppipe=$(cat tmp.tsv | tr ',' '.' | grep -w "${virus}" | grep -w "cooppipe-$virus-consensus.fa" > "${virus}-cooppipe")
  best_tool=$(cat tmp.tsv | tail -n +2 | tr ',' '.'  | grep -w "${virus}" | grep -w --invert-match "cooppipe-$virus-consensus.fa" | head -1 > "${virus}-a")
  
  if [ -s ${virus}-cooppipe ];
    then
    nr_virus=$(echo $nr_virus + 1 |bc -l)
  fi

  
  cd ..
  
  cd NCSD
  cp ../Results/total_stats.tsv .
 
  sort -t$'\t' -k6 -n total_stats.tsv > tmp.tsv
  cooppipe=$(cat tmp.tsv | tr ',' '.' | grep -w "${virus}" | grep -w "cooppipe-$virus-consensus.fa" > "${virus}-cooppipe")
  best_tool=$(cat tmp.tsv | tail -n +2 | tr ',' '.'  | grep -w "${virus}" | grep -w --invert-match "cooppipe-$virus-consensus.fa" | head -1 > "${virus}-a")
  
  cd ..
  
  cd NRC
  cp ../Results/total_stats.tsv .
 
    sort -t$'\t' -k7 -g total_stats.tsv > tmp.tsv
  cooppipe=$(cat tmp.tsv | tr ',' '.' | grep -w "${virus}" | grep -w "cooppipe-$virus-consensus.fa" > "${virus}-cooppipe")
  best_tool=$(cat tmp.tsv | tail -n +2 | tr ',' '.'  | grep -w "${virus}" | grep -w --invert-match "cooppipe-$virus-consensus.fa" | head -1 > "${virus}-a")
 
  cd ..  
done
#
cd avg_identity
find . -type f -empty -delete
rm *.tsv
cd ..
#
cd NCSD
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
list_ncsd=($(ls NCSD))
list_nrc=($(ls  NRC))
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
    ymax = 101
    ymin = 0
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
    aux = 0
    do for [ file in "${list_avg[@]}"]{  
      print file
    
      pos = strstrt(file,"-")
      virus = file[0:pos-1]
     
      set xtics (virus count)
      
      ymax = ymax - offset
      
      if(aux == 0){
        aux = aux + 1;
        plot file using (count):5:13 with labels point pt 7 lc "black" offset char 6,-0.1 notitle
       
      } else {
        plot file using (count):5:13 with labels point pt 3 lc rgb "#009e73" offset char 6,-0.1 notitle
        count = count + spacing_x
        aux = 0
      }
    }
EOF
#
cp *.pdf ../Graphs
cd ../NCSD
#
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Verdade,9'
    set output "NCSD.pdf"
    set datafile separator "\t"
    
    spacing_x = 30
    ymax = 1
    ymin = 0
    offset = ( ymax - ymin ) / 15.0    
    set yrange [ymin:ymax]
    set xrange [0:$nr_virus * spacing_x]
   
    set ytics auto
    set ylabel "NCSD"
    set xlabel "Viruses"
    set multiplot layout 1,1
    set rmargin 5
    set key at screen 1, graph 1 
    
    count = 10
    aux = 0
    do for [ file in "${list_ncsd[@]}"]{  
    
      pos = strstrt(file,"-")
      virus = file[0:pos-1]
     
      set xtics (virus count)
      
      ymax = ymax - offset
      
      if(aux == 0){
        aux = aux + 1;
        plot file using (count):6:13 with labels point pt 7 lc "black" offset char 6,-0.1 notitle
        
       
      } else {
  
        plot file using (count):6:13 with labels point pt 3 lc rgb "#009e73" offset char 6,-0.1 notitle
        count = count + spacing_x
        aux = 0
      }
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
    ymax = 1
    ymin = 0
    offset = ( ymax - ymin ) / 15.0    
    set yrange [ymin:ymax]
    set xrange [0:$nr_virus * spacing_x]
   
    set ytics auto
    set ylabel "NRC"
    set xlabel "Viruses"
    set multiplot layout 1,1
    set rmargin 5
    set key at screen 1, graph 1 
    
    count = 10
    aux = 0
    do for [ file in "${list_nrc[@]}"]{  
    
      pos = strstrt(file,"-")
      virus = file[0:pos-1]
     
      set xtics (virus count)
      
      ymax = ymax - offset
      
      if(aux == 0){
        aux = aux + 1;
        plot file using (count):7:13 with labels point pt 7 lc "black" offset char 6,-0.1 notitle
        
      } else {
        plot file using (count):7:13 with labels point pt 3 lc rgb "#009e73" offset char 6,-0.1 notitle
        count = count + spacing_x
        aux = 0
      }
    }
EOF
#
cp *.pdf ../Graphs
cd ..
#
./Plots_k.sh
#
printf "Finished generating graphs!\n"
#
#set multiplot
#plot file u 13:5 title file with linespoints linestyle count
#rm -rf total_stats.tsv
#, "" u 13:5:5 with labels offset char 0,1
#plot "total_stats.tsv" using 13:6:(sprintf("(%d, %d)", $1, $2)) with labels notitle
#set key at 100., 100.
