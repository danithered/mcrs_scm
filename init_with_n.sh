#!/bin/bash

fromdir=/home/danielred/data/programs/mcrs_chromosome/OUT/A7.1_24_bubble_sok/SAVE/
tempin=IN/tempnbubis

mkdir -p $tempin
no_ok=0

# iterating tru generations
for g in {0..10}
do
	# get list of files from 1st gen
	files=$(ls $fromdir | grep bubble_t${g}_)

	# copy files
	mkdir -p $tempin/$g
	# iterating tru files
	for f in ${files[@]}
	do
		cp ${fromdir}$f $tempin/$g/
		mv $tempin/$g/$f $tempin/$g/bubble_t0_${f}.tsv	
	done

	./mcrscm --par_maxtime 20000 --par_poolsize 1000 --par_splitfrom 50 --par_MN 10 --par_output_interval 10 --par_save_interval 1000 --par_ID A7fromn.3_${f} --par_str_pool IN/str/mappingA7.txt --par_quit 8 --par_bubbles $tempin/$g --par_substitution 0.005 >> fromn.txt

	# if it has died 1
	# if it has survived 0
	# if it was a special conditon: 2

	if [ $? -e 2 ]
	then
		((no_ok++))
		# exit condition
		if [ $no_ok -ge 100 ]
		then
			break
		fi
	fi

	# exit condition
	if [ $no_ok -ge 100 ]
	then
		break
	fi
done

echo Number of simulations with ok condition: $no_ok 
