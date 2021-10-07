#!/bin/bash


# This code will launch a fixed number of replicates for various network sizes

# Note 1: I've put my own path to Simul_prog in line 33, you'll have to put yours. Same for other paths in general.

# Note 2 : About "parallel" (used below, line 36)
# I use here "parallel" to do multiple replicates at the same time, using multicores computing
# "parallel" works like this : you create a text file in which each line is a command to execute.
# You can tell "parallel" how many cores it can use with the -j parameter.
# Then, "parallel" will launch as many tasks as it can, according to the number of cores you allowed.
# You have to install parallel for this to work, and disable the add that appears with the "will cite"

PARAM_TEMPLATE=./template.par
SIMUL_PROG=$HOME/Work/simevolv/bin/Release/Simul_Prog
PARAMFILE_NAME=param.par
rep=20  # Number of replicates for each network size

#-----------------------------8< To be updated
VARIABLE_NAME=GENET_MUTSD 
VARIABLE_NAME2=OUT_CANAL_MUTSD
SHORT_NAME=mutsd1
values=(
    1.0	
    0.1
    0.01
    0.001
    0.0001
    0.00001
)
#->8-------------------------------------------


LAUNCH_FILE=${SHORT_NAME}-launch.sh


if [ -e $LAUNCH_FILE ]
then
	rm -f $LAUNCH_FILE
fi

for j in ${values[@]} # For each mutsd
do
	mydir=${SHORT_NAME}_$j
	
	if [ ! -d $mydir ] 
	then
		mkdir $mydir
	fi
	
	# Create the right parameter file inside each directory
	cat $PARAM_TEMPLATE | sed "s/${VARIABLE_NAME}.*/${VARIABLE_NAME}\t$j/" | sed "s/${VARIABLE_NAME2}.*/${VARIABLE_NAME2}\t$j/" > $mydir/$PARAMFILE_NAME 

	for i in `seq 1 $rep`; # For each replicate
	do
		echo $SIMUL_PROG -p $mydir/$PARAMFILE_NAME -o $mydir/simulCov${i}.txt >> $LAUNCH_FILE
		# Writing this command in a text file
	done 
	
	#~ parallel -a ./launcher.sh -j $par # Launching parallel
	
done
	
# This can be done in a better / faster way, but this way seems better to me for learning at first
# Have fun ! :)
