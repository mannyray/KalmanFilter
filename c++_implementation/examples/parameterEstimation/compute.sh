#!/bin/bash

echo "compiling"
g++ basicExampleNonlinear.cpp -o script

actualProcessNoise=0.1;
actualSensorNoise=0.01;
for assumedProcessNoise in $(seq 0.05 0.01 0.25); do
	echo "ASSUMED PROCESS NOISE: " $assumedProcessNoise
	#for assumedSensorNoise in $(seq 0.05 0.01 0.1); do
		#echo "assumed sensor noise: " $assumedSensorNoise
		fileSave="assumedP_"$assumedProcessNoise"assumedR"$actualSensorNoise	
		./script $actualProcessNoise $actualSensorNoise $assumedProcessNoise $actualSensorNoise $fileSave
	#done
done




