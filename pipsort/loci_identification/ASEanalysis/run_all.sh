infile=inputs.csv
infile=remaining.csv
while IFS=$', ' read -r -a myArray
do
 echo ${myArray[0]}
 echo ${myArray[1]}
 echo ${myArray[2]}
 echo ${myArray[3]}
 echo ${myArray[4]}
 echo ${myArray[5]}
 #bash full_local_ancestry.sh ${myArray[0]} ${myArray[1]} ${myArray[2]} ${myArray[3]} ${myArray[4]} ${myArray[5]} ${myArray[6]}
 #bash interaction_gwas.sh ${myArray[0]} ${myArray[1]} ${myArray[2]} ${myArray[3]} ${myArray[4]} ${myArray[5]} ${myArray[6]}
 bash no_heterozygous_interaction.sh ${myArray[0]} ${myArray[1]} ${myArray[2]} ${myArray[3]} ${myArray[4]} ${myArray[5]} ${myArray[6]}
done < $infile
