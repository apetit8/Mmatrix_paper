#!/bin/bash


SIMUL_PROG="/shared/projects/evoplanet/Software/simevolv/bin/Release/Simul_Prog"
rep=30  # Number of replicates for each network size
paramfiles=`ls simul/fig_4/*/*.par`
SHORT_NAME=All_fig
LAUNCH_FILE=${SHORT_NAME}-launch.sh
PARAMFILE_NAME=param.par


if [ -e $LAUNCH_FILE ]
then
	rm -f $LAUNCH_FILE
fi


for j in $paramfiles
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