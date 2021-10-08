#!/bin/bash 

# read in template file
infile=$1

n=1
pos=()
base=()
while read line
do
	if [ $n -eq 1 ]
	then
		str=$line
	else
		pos+=($(echo $line | awk '{print $1}'))
		base+=($(echo $line | awk '{print $2}'))
	fi
	n=$((n+1))
done < $infile

#echo ${base[2]}
#grep "(.(......))" OUT/randseqs.txt | grep "UUAGAACGACCA" | awk '{print S2}'

# read in sequences and structures
vals=($( grep -F $str OUT/randseqs.txt ))

# test for subrules
nook=0
for ((i = 0 ; i < ${#vals[@]} ; i += 4 ))
do
	#calculating where does the structure starts
	rest=${vals[$((i+1))]#*$str}
	from=$(( ${#vals[$i]} - ${#rest} - ${#str} ))
	restseq=${vals[$i]:$from}
	
	passed=1 
	
	#testing subrules
	for ((b = 0 ; b < ${#base[@]} ; b++))
	do
		if [ $(expr substr $restseq ${pos[$b]} 1) != ${base[$b]} ]
		then 
			#echo testing ${vals[$i]} in pos ${pos[$b]} for base ${base[$b]} 
		#else
			passed=0
			break
		fi
	done
	
	if [ $passed == 1 ]
	then
		#echo tested ${vals[$i]} 
		#echo ${vals[$((i+1))]:$from}
		#echo $restseq
		nook=$((nook + 1))
	fi
	
done

#resutls
echo str found in $(( ${#vals[@]} / 4)) sequences and subrules applied for $nook
