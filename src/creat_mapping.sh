#!/bin/bash 

no_strs=$1

#first line of mapping file (number of activities)
echo $no_strs

#going thru directories
for(( d=1 ; d <= $no_strs ; d++ ))
do	
	#find str files in directory
	strfiles=($(find IN/str/${d}/str*))
	
	#looping thru those files
	for s in ${strfiles[@]}
	do
		startlines=($(awk '/\./ {print FNR}' $s) $(( $(wc -l < $s) + 1 )) )
		
		for n in $(seq 0 $(( ${#startlines[@]} -2 )))
		do
			#print rule
			printf "1 "
			sed -n ${startlines[n]}p $s
			
			#print subrules
			awk -v s="${startlines[n]}" -v e="${startlines[n+1]}" 'NR>s&&NR<e' $s | sed ':a;N;$!ba;s/\n/ /g'
			
			#print activities for subrules
			for (( a = 1 ; a <= $no_strs ; a++ ))
			do
				if [ $a = $d ]
				then
					printf "%d " 1
				else
					printf "%d " 0
				fi
			done
			echo
		done
	done
done
