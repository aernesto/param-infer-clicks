#!/bin/bash
set -e
# recall produce_figs(dbname, dsetname, numtrials, numsamples,init_sample, last_sample, savelocation)
filename="data/S6lr30h20T2tr1000sp5000.h5"
dset="/lr30hr97.84h20T2"
ntrials=1000
nsamp=5000
init_h=0
last_h=40
init_g=40
last_g=80
save="figs/progressionS6lr30h20T2/"
matlab -nodisplay -r "produce_figs('$filename', '$dset', $ntrials, $init_g, $last_g, $nsamp, $init_h, $last_h, '$save')" &> logs/produce_figs_log.txt&
