#!/bin/bash


SIMUL_PROG="/shared/projects/evoplanet/Software/simevolv/bin/Release/Simul_Prog"
rep=1  # Number of replicates for each network size
templates=`ls simul/fig_6cor/*/*.par`
SHORT_NAME=fig_1abc
LAUNCH_FILE=${SHORT_NAME}-launch.sh
PARAMFILE_NAME=param.par


if [ -e $LAUNCH_FILE ]
then
	rm -f $LAUNCH_FILE
fi


for j in $templates
do
	mydir=${j%%.par*}
	
	if [ ! -d $mydir ] 
	then
		mkdir $mydir # Create a directory for this paramfile
	fi
	
	# Create the right parameter file inside each directory
	cat $j > $mydir/$PARAMFILE_NAME 
	rm $j

	for i in `seq 1 $rep`; # For each replicate
	do
		echo $SIMUL_PROG -p $mydir/$PARAMFILE_NAME -o $mydir/simulAngle${i}.txt >> $LAUNCH_FILE
		# Writing this command in a text file
	done 
	
done
