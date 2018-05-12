#!/bin/bash
set -e
# this script does the following:
# 1/ creates a hdf5 file with Python 3
# 2/ computes the decision data from the linear model with Python 3
# 3/ computes the decision data from the nonlinear model with MATLAB
# 4/ produces figures of proportion of compatible trials with MATLAB

# SET PARAMETERS
S=3;
lowrate=1;
highrate="`python3 get_lambda_high.py $lowrate $S`";
kappa="`python3 get_kappa.py $lowrate $highrate`";
file_substr="ttt";
filename="data/$file_substr.h5";
logfile="logs/create_$file_substr""_log.txt"
h=1;
T=2;
groupname="/lr$lowrate""hr$highrate""h$h""T$T/";
ntrials=100;
nsamples=50;

# Uncomment following block for stdout checks
: <<'END'
echo
echo "S: "$S 
echo "lr: "$lowrate 
echo "hr: "$highrate 
echo "kappa: "$kappa 
echo "filename: "$filename 
echo "h: "$h 
echo "T: "$T 
echo "groupname: "$groupname
echo
END

echo "STEP 1-2" > $logfile
python3 create_db.py "$lowrate" "$S" "$h" "$T" "$filename" "$ntrials" "$nsamples" &>> $logfile
echo "STEP 3" >> $logfile
matlab -nodesktop -nodisplay -r "fill_nonlin_dec($lowrate,$highrate,$kappa,$h,$T,'$filename',$ntrials,$nsamples)" &>> $logfile &
echo "STEP 4" >> $logfile
#matlab -nodesktop -nodisplay < produce_figs.m &>> $logfile
