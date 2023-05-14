#!/bin/bash

fromdir=/home/danielred/data/programs/mcrs_chromosome/OUT/A7.1_24_bubble_sok/SAVE/

mkdir -p IN/tempbubis
no_ok=0

# iterating tru generations
for g in {0..10}
do
	# get list of files from 1st gen
	files=$(ls $fromdir | grep bubble_t${g}_)
	

	# iterating tru files
	for f in ${files[@]}
	do
		# copy file
		mkdir IN/tempbubis/$f
		cp ${fromdir}$f IN/tempbubis/$f/
		mv IN/tempbubis/$f/$f IN/tempbubis/$f/bubble_t0_1.tsv	

		# running simulation
		./mcrscm --par_maxtime 20000 --par_poolsize 1000 --par_splitfrom 50 --par_MN 10 --par_output_interval 10 --par_save_interval 1000 --par_ID A7fromone.3_${f} --par_str_pool IN/str/mappingA7.txt --par_quit 4 --par_bubbles IN/tempbubis/$f --par_substitution 0.005 >> fromone.txt

		# if it has died 1
		# if it has survived 0
		# if it was a special conditon: 2

		if [ $? -e 2 ] # 2: special condition: split
		then
			((no_ok++))
			# exit condition
			if [ $no_ok -ge 100 ]
			then
				break
			fi
		fi
	done

	# exit condition
	if [ $no_ok -ge 100 ]
	then
		break
	fi
done

echo Number of simulations with ok condition: $no_ok 
