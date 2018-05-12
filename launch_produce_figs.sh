#!/bin/bash
set -e
# recall produce_figs(dbname, dsetname, numtrials, numsamples,init_sample, last_sample, savelocation)
filename="data/S3lr2h1T2tr5Ksp10K.h5"
dset="/lr2hr14h1T2"
ntrials=5000
nsamp=10000
initsamp=0
lastsamp=40
save="figs/"
matlab -nodisplay -r "produce_figs('$filename', '$dset', $ntrials, $nsamp,$initsamp, $lastsamp, '$save')" &> logs/produce_figs_log.txt&
