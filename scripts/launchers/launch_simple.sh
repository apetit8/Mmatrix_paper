#!/bin/bash


# This code will launch a fixed number of replicates for various network sizes

# Note 1: I've put my own path to Simul_prog in line 33, you'll have to put yours. Same for other paths in general.

# Note 2 : About "parallel" (used below, line 36)
# I use here "parallel" to do multiple replicates at the same time, using multicores computing
# "parallel" works like this : you create a text file in which each line is a command to execute.
# You can tell "parallel" how many cores it can use with the -j parameter.
# Then, "parallel" will launch as many tasks as it can, according to the number of cores you allowed.
# You have to install parallel for this to work, and disable the add that appears with the "will cite"





rep=100  # Number of replicates for each network size
par=24   # Number of cores to use




mkdir $j"simu" # Create a directory
cd $j"simu"
cat ../param*.par | sed "s/GENET_NBLOC/GENET_NBLOC $j/g" > param.par # Create the right parameter file inside each directory

	for i in `seq 1 $rep`; # For each replicate
	do
		echo ../../../../simevolv/bin/Release/Simul_Prog -p ./param.par -o simulCov${i}.txt # Writing this command in a text file
	done > ./launcher.sh
	
	parallel -a ./launcher.sh -j $par # Launching parallel


cd ..
	

	
# This can be done in a better / faster way, but this way seems better to me for learning at first
# Have fun ! :)
