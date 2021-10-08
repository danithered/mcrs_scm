#!/bin/bash 

sleep 1

#saving arguments
if [ "$#" -lt 5 ]; then
    echo "Illegal number of parameters"
fi

no_strs=$1
ncol=$2
nrow=$3
dens=$4
nc=$5

#calculate stuff
length=$((ncol*nrow))

#going thru directories to get number of all sequences
Nall=0
for(( d=1 ; d <= $no_strs ; d++ ))
do	
	#find seqfiles in directory
	file=($(find IN/str/$d/randseqs* | head -n 1))
	Nall=$(( Nall + $(wc -l < $file) ))
	#echo Number of all strs: $Nall
done
#enrichment=$(( length*dens/Nall/100)) 
#echo Enrichment is $enrichment

#going thru directories to sample it
for(( d=1 ; d <= $no_strs ; d++ ))
do	
	#find seqfiles in directory
	#file=($(find IN/str/$d/randseqs_ea* | head -n 1))
	file=($(find IN/str/$d/ -name randseqs_ea* -and ! -name *_compl.txt))
	
	if [ "$nc" == "nc" ] 
	then
		# in case there is no complementer generated
		perfile=$(echo "scale=10; x=$length*$dens/$Nall/100 * $(wc -l < $file); scale=0; x/1" | bc )
		shuf -r -n $perfile -r $file >> temp.txt
	else
		# in case of complementers
		compl_file=($(find IN/str/$d/randseqs_ea*_compl.txt | head -n 1))
		perfile=$(echo "scale=10; x=$length*$dens/$Nall/100 * $(wc -l < $file); scale=0; x/2" | bc )
		shuf -r -n $perfile -r $file >> temp.txt
		shuf -r -n $perfile -r $compl_file >> temp.txt
	fi
done

no_lines=$( wc -l < temp.txt )
no_empty=$((length-no_lines))
#echo adding $no_empty empty lines
for (( i=no_empty; i > 0; i-- ))
do
	echo "N N 0" >> temp.txt
done

shuf temp.txt
rm temp.txt
