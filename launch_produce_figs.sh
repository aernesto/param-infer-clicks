#!/bin/bash
set -e
# recall produce_figs(dbname, dsetname, numtrials, numsamples,init_sample, last_sample, savelocation)
filename="data/S3lr2h1T2tr500sp5K.h5"
dset="/lr2hr14h1T2"
ntrials=500
nsamp=5000
initsamp=0
lastsamp=40
save="figs/progressionS3lr2h1T2/"
matlab -nodisplay -r "produce_figs('$filename', '$dset', $ntrials, $nsamp,$initsamp, $lastsamp, '$save')" &> logs/produce_figs_log.txt&
