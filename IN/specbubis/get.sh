#!/bin/bash

files=$(awk '{print $1;}' descriptions.txt)

count=1

for file in ${files[@]}
do
	# prepare targetdir
	rm -rf ${count}
	mkdir ${count}
	
	cp $file ${count}/
	count=$((${count}+1))
done
