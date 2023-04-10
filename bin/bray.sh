#!/bin/bash
#
#Usage: ./bray.sh otuTable.csv
#
# This script is inteded to perform a bary-curtis dissimilarity matrix calculation between all pairs of metagenomes. When the matrix is ready the script
#will automatically use the matrix to represent the distances with a UPGMA tree. The input table can be derived from the OTU.table, however the table must
#be comma (,) delimited and have a header consisting of an empty first cell (deleting first column name) followed by all the samples names.
#
# For tree building it requires neighbor from the phylip package (https://evolution.genetics.washington.edu/phylip.html).
#
#
n=$(head -n1 $1 | tr -s ',' '\n' | wc -l | tr -d ' ')
for i in $(head -n1 $1 | tr -s ',' '\n' | tr -d '"' | sed 's/\.kaiju\.tsv//'); do echo $i | sed 's/'"$i"'/'"$i"'         /' | cut -c -9 ; done | sed '1s/^/\'$'\n/' > 1names.tmp
for i in $(seq 2 $((n-1))); do 
    printf '\n%.0s' $(seq $i) > "$i".tmp
    for j in $(seq $((i+1)) $n); do
        #echo $i y $j
        met1=$(awk -F ',' '{print $'$i'}' $1 | head -n1 | sed 's/\.kaiju\.tsv//')
        met2=$(awk -F ',' '{print $'$j'}' $1 | head -n1 | sed 's/\.kaiju\.tsv//')  
        #echo "$met1 y $met2" 
        min_sum=$(awk -F ',' '{if ($'$i' < $'$j') {print $'$i'} else {print $'$j'} }' $1 | awk -F ',' '{s+=$1}END{print s}')
        #echo "minimal sum is $min_sum"
        times2=$(( 2*min_sum ))
        #echo "minimal sum times 2 is $times2"
        sum1=$(awk  -F ',' '{s+=$'$i'}END{print s}' $1)
        #echo "Sum1 is $sum1"
        sum2=$(awk -F ',' '{s+=$'$j'}END{print s}' $1)
        #echo "Sum2 is $sum2"
        total_sum=$(( sum1+sum2 ))
        #echo "total sum is $total_sum"
        division=$(echo $times2 / $total_sum | bc -l)
        #echo "division is $division"
        distance=$(echo 1 - $division | bc)
        echo $distance
    done >> "$i".tmp
done
n2=$(( n-1 ))
ls *.tmp | sort -n | xargs paste | sed '1d' | sed '1s/^/'"$n2"'\'$'\n/' | sed $'s/\t/ /g' > "${1%Table.csv}_matrix.txt"
rm *.tmp

echo ""${1%Table.csv}_matrix.txt"
N
L
Y
" > ${1%Table.csv}.input
neighbor < ${1%Table.csv}.input > "${1%Table.csv}_screen.out"
mv outfile "${1%Table.csv}_matrix.out"
mv outtree "${1%Table.csv}_matrix.tree"
