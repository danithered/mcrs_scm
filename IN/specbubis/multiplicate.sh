#!/bin/bash

for i in */
do 
	files=$i/*
	count=1
	
	# copy files to tempdir
	mkdir -p tempdir
	for f in ${files[@]}
	do
		cp $f tempdir/$(basename $f)
	done	
	
	# copy files a number of times
	while [ $count -le $1 ]
	do
		# re-use files
		for f in tempdir/*
		do
			if [ $count -gt $1 ]
			then
				break
			fi
			cp $f $i/bubble_t0_${count}.tsv
			count=$((${count}+1))
		done
	done
	
	rm -r tempdir
done


