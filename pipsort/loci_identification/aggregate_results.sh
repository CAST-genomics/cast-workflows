logfile=$1
infile=regions.txt
while IFS=$' ' read -r -a myArray
do
 echo ${myArray[0]}
 echo ${myArray[1]}
 echo ${myArray[2]}
 cat ${myArray[0]}_${myArray[1]}_${myArray[2]}_cluster/${myArray[0]}_${myArray[1]}_${myArray[2]}_cluster_costs >> $logfile
done < $infile
