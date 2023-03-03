#!/bin/bash

cd OUT
odirs=( $(ls -d -1 */) )
cd ..

for od in ${odirs[@]}
do
	FILE=OUT/${od}output.csv
	
	#does output exists?
	if test -f "$FILE"
	then
		IFS=';' read -r -a outp <<< "$(tail -n 1 $FILE)" # read in last line of output
		if [ ${outp[1]} -gt 100 ]
		then
			echo -e $od '\t'is alive \(${outp[1]}\) at time ${outp[0]} >> isalive.temp
		#else
			#echo -e $od '\t'is dead >> isalive.temp
		fi
	fi
done

column -t -s' ' isalive.temp
rm isalive.temp
