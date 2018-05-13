#!/bin/bash
set -e
# recall produce_figs(dbname, dsetname, numtrials, numsamples,init_sample, last_sample, savelocation)
filename="data/S4lr2h4T2tr5Ksp5K.h5"
dset="/lr2hr21.31h4T2"
ntrials=5000
nsamp=5000
init_h=0
last_h=40
init_g=320
last_g=360
save="figs/progressionS6lr30h20T2/"
matlab -nodisplay -r "produce_figs('$filename', '$dset', $ntrials, $init_g, $last_g, $nsamp, $init_h, $last_h, '$save')" &> logs/produce_figs_log.txt&
