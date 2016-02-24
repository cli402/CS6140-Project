#!/bin/bash


graphFiles=`ls ./data/ | grep .gr`

for graph in $graphFiles
do	 	
	filename=`echo $graph | cut -d'.' -f1`
	echo $graph $filename
	python ./src/run_experiments.py $graph $filename.extra ./results/$filename_output.txt

done
