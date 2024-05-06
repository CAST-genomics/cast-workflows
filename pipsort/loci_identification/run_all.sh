infile=platelet_count_1MB_loci.txt
while IFS=$' ' read -r -a myArray
do
 echo ${myArray[0]}
 echo ${myArray[1]}
 echo ${myArray[2]}
 echo ${myArray[3]}
 bash run_pipsort.sh ${myArray[0]} ${myArray[1]} ${myArray[2]} ${myArray[3]} EUR_WHITE.csv AFR_BLACK.csv platelet_count
done < $infile
