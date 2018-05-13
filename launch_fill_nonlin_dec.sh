#!/bin/bash
set -e
# recall fill_nonlin_dec(lr, hr, kappa, h, td, dbname, ntrials, nsamples)
#lowrate="30"
lowrate="2"
highrate="21.313708498984759"
#highrate="97.8397431775"
kap="2.3662032766408734"
#kap="1.1821334847899723"
#haz="20"
haz="4"
T="2"
filename="data/S4lr2h4T2tr5Ksp5K.h5"
#filename="data/S6lr30h20T2tr500sp1000.h5"
#dset="/lr30hr97.84h20T2"
ntrials=5000
nsamp=5000
matlab -nodisplay -r "fill_nonlin_dec($lowrate,$highrate,$kap,$haz,$T,'$filename',$ntrials,$nsamp)" &> logs/fill_nonlin_dec_log_S4lr2h4.txt&
