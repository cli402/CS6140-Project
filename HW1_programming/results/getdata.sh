#!/bin/bash

file=$(ls r*)
for f in $file
 do 	
	total=0
	cat $f | awk '{total = total + $2}END{print total '\n'}' >>time.txt
done
