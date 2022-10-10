#!/bin/bash 

sleep 1

#saving arguments
no_strs=$1
ncol=$2
nrow=$3
dens=$4

#calculate stuff
length=$((ncol*nrow))
perfile=$((length*dens/100/no_strs))

#echo $perfile

#going thru directories
for(( d=1 ; d <= $no_strs ; d++ ))
do	
	#find seqfiles in directory
	#file=($(find IN/str/$d/randseqs* | head -n 1))
	file=($(find IN/str/$d/ -name randseqs_ea* -and ! -name *_compl.txt))


	#echo $file $perfile
	#pwd
	shuf -n $perfile -r $file >> temp.txt
done

no_lines=$( wc -l < temp.txt )
no_empty=$((length-no_lines))
for (( i=no_empty; i > 0; i-- ))
do
	echo "N N 0" >> temp.txt
done

shuf temp.txt
rm temp.txt
