#!/bin/bash


SIMUL_PROG=$HOME/Work/simevolv/bin/Release/Simul_Prog
rep=20  # Number of replicates for each network size
#par=24   # Number of cores to use

templates=`ls ./*.par`
SHORT_NAME=ecc
LAUNCH_FILE=${SHORT_NAME}-launch.sh
PARAMFILE_NAME=param.par


if [ -e $LAUNCH_FILE ]
then
	rm -f $LAUNCH_FILE
fi


for j in $templates
do
	short=${j##*ecc}
	mydir=${SHORT_NAME}_$short
	
	if [ ! -d $mydir ] 
	then
		mkdir $mydir # Create a directory for this size of network
	fi
	
	# Create the right parameter file inside each directory
	cat $j > $mydir/$PARAMFILE_NAME 

	for i in `seq 1 $rep`; # For each replicate
	do
		echo $SIMUL_PROG -p $mydir/$PARAMFILE_NAME -o $mydir/simulEcc${i}.txt >> $LAUNCH_FILE
		# Writing this command in a text file
	done 
	
	#~ parallel -a ./launcher.sh -j $par # Launching parallel
done
