#!/bin/bash

datafile=`ls ./data | grep .txt`

for data in ${datafile[*]}
do 
	#echo $data	
	python ./hw3.py $data
done
