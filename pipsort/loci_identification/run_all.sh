#infile=platelet_count_1MB_loci.txt
#infile=remaining_loci.txt
#infile=platelet_count_loci_1MB_AFR_NOT_AFR.txt
#infile=small_loci.txt
infile=white_blood_cell_count_NOT_AFR_AFR_1MB_loci.txt
#infile=egfr_ckdepi_1MB_loci.txt
while IFS=$' ' read -r -a myArray
do
 echo ${myArray[0]}
 echo ${myArray[1]}
 echo ${myArray[2]}
 echo ${myArray[3]}
 #bash run_pipsort.sh ${myArray[0]} ${myArray[1]} ${myArray[2]} ${myArray[3]} EUR_WHITE.csv AFR_BLACK.csv platelet_count
 #bash run_pipsort.sh ${myArray[0]} ${myArray[1]} ${myArray[2]} ${myArray[3]} NOT_AFR_BLACK.csv AFR_BLACK.csv platelet_count
 #bash run_pipsort.sh ${myArray[0]} ${myArray[1]} ${myArray[2]} ${myArray[3]} NOT_AFR_BLACK.csv AFR_BLACK.csv egfr_ckdepi
 bash run_pipsort.sh ${myArray[0]} ${myArray[1]} ${myArray[2]} ${myArray[3]} NOT_AFR_BLACK.csv AFR_BLACK.csv white_blood_cell_count
done < $infile
