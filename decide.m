ntrials=100000;
ncols=2;
ncols_write=10001;
all_h = [1, linspace(0,40,ncols_write-1)];
td = 2;

% S2 data
dbname = '/storage/adrian/data_S_2_5.h5';
dsetname = '/lr1hr6.46h1T2/trials';
dsetname_decision = '/lr1hr6.46h1T2/decision_nonlin';
trial_data = h5read(dbname,dsetname,[1 1],[ncols ntrials]);
kappa = log(6.4641016151377544);
for i = 1:ntrials
    ls = trial_data{1,i};
    rs = trial_data{2,i};
    dec = zeros(ncols_write,1);
    parfor tn=1:ncols_write
        dec(tn) = decide1(td, ls, rs, kappa, all_h(tn), 0);
    end
    start=[1 i];
    count=[ncols_write 1];
    h5write(dbname, dsetname_decision, dec, start, count)
end
