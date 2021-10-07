PARAM_TEMPLATE=./test_corr_1_2.par
SIMUL_PROG=$HOME/simevolv/bin/Release/Simul_Prog
PARAMFILE_NAME=param.par
rep=20  # Number of replicates for each network size

#-----------------------------8< To be updated
VARIABLE_NAME=FITNESS_CORRELATION
SHORT_NAME=test
values=(
    "-0.8 0 0 0 0 0"
    "0.8 0 0 0 0 0"
    "0 0 0 0 0 0"
)
#->8-------------------------------------------


LAUNCH_FILE=${SHORT_NAME}-launch.sh


if [ -e $LAUNCH_FILE ]
then
	rm -f $LAUNCH_FILE
fi

for j in ${values[@]}
do
	mydir=${SHORT_NAME}_$j
	
	if [ ! -d $mydir ] 
	then
		mkdir $mydir # Create a directory for this size of network
	fi
	
	# Create the right parameter file inside each directory
	cat $PARAM_TEMPLATE | sed "s/${VARIABLE_NAME}/${VARIABLE_NAME}\t$j/" > $mydir/$PARAMFILE_NAME 

	for i in `seq 1 $rep`; # For each replicate
	do
		echo $SIMUL_PROG -p $mydir/$PARAMFILE_NAME -o $mydir/simulCov${i}.txt >> $LAUNCH_FILE
		# Writing this command in a text file
	done 
	
	#~ parallel -a ./launcher.sh -j $par # Launching parallel
	
done
