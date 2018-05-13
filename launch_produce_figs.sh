#!/bin/bash
set -e
# recall produce_figs(dbname, dsetname, numtrials, numsamples,init_sample, last_sample, savelocation)
filename="data/S4lr2h4T2tr5Ksp5K.h5"
dset="/lr2hr21.31h4T2"
ntrials=5000
nsamp=5000
initsamp=0
lastsamp=40
save="figs/progressionS4lr2h4T2/"
matlab -nodisplay -r "produce_figs('$filename', '$dset', $ntrials, $nsamp,$initsamp, $lastsamp, '$save')" &> logs/produce_figs_log.txt&
