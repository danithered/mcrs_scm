#!/bin/bash
for file in ./OUT/randseqs_ea* 
do
	#echo $file
	head $file -n 500
done 
