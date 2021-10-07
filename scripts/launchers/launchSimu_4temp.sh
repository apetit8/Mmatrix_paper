PARAM_TEMPLATE=./templatec0.5_0.5_10_3.par
SIMUL_PROG=$HOME/Work/simevolv/bin/Release/Simul_Prog
PARAMFILE_NAME=param.par
rep=20  # Number of replicates for each network size

#-----------------------------8< To be updated
VARIABLE_NAME=FITNESS_CORRELATION
SHORT_NAME=c0.5_opt0.5_s10_3_corr
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

for j in ${values[@]} # For each size of network from 2 to maxSize
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
	
# This can be done in a better / faster way, but this way seems better to me for learning at first
# Have fun ! :)
########################################################################

PARAM_TEMPLATE2=./templatec0.5_0.8_10_3.par
rep=20  # Number of replicates for each network size


#-----------------------------8< To be updated
VARIABLE_NAME=FITNESS_CORRELATION
SHORT_NAME2=c0.5_opt0.8_s10_3_corr
values2=(
    "-0.8 0 0 0 0 0"
    "0.8 0 0 0 0 0"
    "0 0 0 0 0 0"
)
#->8-------------------------------------------


LAUNCH_FILE=${SHORT_NAME2}-launch.sh


if [ -e $LAUNCH_FILE ]
then
	rm -f $LAUNCH_FILE
fi

for j in ${values2[@]} # For each size of network from 2 to maxSize
do
	mydir=${SHORT_NAME2}_$j
	
	if [ ! -d $mydir ] 
	then
		mkdir $mydir # Create a directory for this size of network
	fi
	
	# Create the right parameter file inside each directory
	cat $PARAM_TEMPLATE2 | sed "s/${VARIABLE_NAME}/${VARIABLE_NAME}\t$j/" > $mydir/$PARAMFILE_NAME 

	for i in `seq 1 $rep`; # For each replicate
	do
		echo $SIMUL_PROG -p $mydir/$PARAMFILE_NAME -o $mydir/simulCov${i}.txt >> $LAUNCH_FILE
		# Writing this command in a text file
	done 
	
	#~ parallel -a ./launcher.sh -j $par # Launching parallel
	
done

#########################################################################

PARAM_TEMPLATE3=./template0.8_10_6.par
rep=20  # Number of replicates for each network size


#-----------------------------8< To be updated
VARIABLE_NAME=FITNESS_CORRELATION
SHORT_NAME3=c0.5_opt0.8_s10_6_corr
values3=(
    "-0.8 0 0 0 0 0"
    "0.8 0 0 0 0 0"
    "0 0 0 0 0 0"
)
#->8-------------------------------------------


LAUNCH_FILE=${SHORT_NAME3}-launch.sh


if [ -e $LAUNCH_FILE ]
then
	rm -f $LAUNCH_FILE
fi

for j in ${values3[@]} # For each size of network from 2 to maxSize
do
	mydir=${SHORT_NAME3}_$j
	
	if [ ! -d $mydir ] 
	then
		mkdir $mydir # Create a directory for this size of network
	fi
	
	# Create the right parameter file inside each directory
	cat $PARAM_TEMPLATE3 | sed "s/${VARIABLE_NAME}/${VARIABLE_NAME}\t$j/" > $mydir/$PARAMFILE_NAME 

	for i in `seq 1 $rep`; # For each replicate
	do
		echo $SIMUL_PROG -p $mydir/$PARAMFILE_NAME -o $mydir/simulCov${i}.txt >> $LAUNCH_FILE
		# Writing this command in a text file
	done 
	
	#~ parallel -a ./launcher.sh -j $par # Launching parallel
	
done

########################################################################

PARAM_TEMPLATE4=./template0.5_10_6.par
rep=20  # Number of replicates for each network size


#-----------------------------8< To be updated
VARIABLE_NAME=FITNESS_CORRELATION
SHORT_NAME4=c0.5_opt0.5_s10_6_corr
values4=(
    "-0.8 0 0 0 0 0"
    "0.8 0 0 0 0 0"
    "0 0 0 0 0 0"
)
#->8-------------------------------------------


LAUNCH_FILE=${SHORT_NAME4}-launch.sh


if [ -e $LAUNCH_FILE ]
then
	rm -f $LAUNCH_FILE
fi

for j in ${values4[@]} # For each size of network from 2 to maxSize
do
	mydir=${SHORT_NAME4}_$j
	
	if [ ! -d $mydir ] 
	then
		mkdir $mydir # Create a directory for this size of network
	fi
	
	# Create the right parameter file inside each directory
	cat $PARAM_TEMPLATE4 | sed "s/${VARIABLE_NAME}/${VARIABLE_NAME}\t$j/" > $mydir/$PARAMFILE_NAME 

	for i in `seq 1 $rep`; # For each replicate
	do
		echo $SIMUL_PROG -p $mydir/$PARAMFILE_NAME -o $mydir/simulCov${i}.txt >> $LAUNCH_FILE
		# Writing this command in a text file
	done 
	
	#~ parallel -a ./launcher.sh -j $par # Launching parallel
	
done
