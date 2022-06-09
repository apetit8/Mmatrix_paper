#!/bin/bash


SIMUL_PROG="/shared/projects/transp_horizon/Software/simevolv/bin/Release/Simul_Prog"
cd ../simul
rep=30  # Number of replicates for each network size
templates=`find -type f -wholename "./*.par"`
SHORT_NAME=fig_
LAUNCH_FILE=${SHORT_NAME}-launch.sh
PARAMFILE_NAME=param.par


if [ -e $LAUNCH_FILE ]
then
	rm -f $LAUNCH_FILE
fi


for j in $templates
do
	short=${j##*simu}
	mydir=$j$short
	
	if [ ! -d $mydir ] 
	then
		mkdir $mydir # Create a directory for this paramfile
	fi
	
	# Create the right parameter file inside each directory
	cat $j > $mydir/$PARAMFILE_NAME 

	for i in `seq 1 $rep`; # For each replicate
	do
		echo $SIMUL_PROG -p $mydir/$PARAMFILE_NAME -o $mydir/simulAngle${i}.txt >> $LAUNCH_FILE
		# Writing this command in a text file
	done 
	
	#~ parallel -a ./launcher.sh -j $par # Launching parallel
done
