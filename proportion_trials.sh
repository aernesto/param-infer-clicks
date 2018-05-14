#!/bin/bash
set -e
# this script does the following:
# 1/ creates a hdf5 file with Python 3
# 2/ computes the decision data from the linear model with Python 3
# 3/ computes the decision data from the nonlinear model with MATLAB
# 4/ produces figures of proportion of compatible trials with MATLAB

# SET PARAMETERS
S=6
lowrate=30
h=20
T=2
ntrials=100000
nsamples=3
gamma_init=0
gamma_end=40
h_init=0
h_end=40
highrate="`python3 get_lambda_high.py $lowrate $S`"
kappa="`python3 get_kappa.py $lowrate $highrate`"
file_substr="S$S""lr$lowrate""h$h""T$T""tr$ntrials""sp$nsamples"
filename="data/$file_substr.h5"
logfile="logs/create_$file_substr""_log.txt"
groupname="`python3 build_group_name.py $lowrate $highrate $h $T`"
saveimfolder="figs/progressionS$S""lr$lowrate""h$h""T$T"
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
python3 create_db.py "$lowrate" "$S" "$h" "$T" "$filename" "$ntrials" "$nsamples" "$gamma_init" "$gamma_end" &>> $logfile
echo "STEP 3 -- AGAIN" >> $logfile
matlab -nodisplay -r "fill_nonlin_dec($lowrate,$highrate,$kappa,$h,$T,'$filename',$ntrials,$nsamples,[$h_init,$h_end])" &>> $logfile &
echo -n "`date`," >> status_db.txt
echo -n "$filename," >> status_db.txt
echo -n "$groupname," >> status_db.txt
echo "hr$highrate" >> status_db.txt
#echo "STEP 4" >> $logfile
# following line doesn't work as such
#matlab -nodisplay -r "produce_figs('$filename', '$groupname', $ntrials,$nsamples,0, 40, '$saveimfolder')" &>> $logfile &
