infile=inputs.csv
infile=remaining.csv
while IFS=$', ' read -r -a myArray
do
 echo ${myArray[0]} #chr
 echo ${myArray[1]} #from
 echo ${myArray[2]} #to
 echo ${myArray[3]} #snpposhg38
 #echo ${myArray[4]} #snpposhg19
 echo ${myArray[5]} #phenname
 echo ${myArray[6]} #phenocovar file
 bash quartile_strat_local_ancestry.sh ${myArray[0]} ${myArray[1]} ${myArray[2]} ${myArray[3]} ${myArray[5]} ${myArray[6]}
done < $infile
