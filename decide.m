%{
	Content of srvr_data_1.h5
	--------------------------
	lr15hr101.258639865h1T2
	lr15hr101.258639865h1T2/decision_lin
	lr15hr101.258639865h1T2/decision_nonlin
	lr15hr101.258639865h1T2/trial_info
	lr15hr101.258639865h1T2/trials
	lr15hr17.8664640249h1T2
	lr15hr17.8664640249h1T2/decision_lin
	lr15hr17.8664640249h1T2/decision_nonlin
	lr15hr17.8664640249h1T2/trial_info
	lr15hr17.8664640249h1T2/trials
	lr15hr36.5367250374h1T2
	lr15hr36.5367250374h1T2/decision_lin
	lr15hr36.5367250374h1T2/decision_nonlin
	lr15hr36.5367250374h1T2/trial_info
	lr15hr36.5367250374h1T2/trials
%}


ntrials=100000;
ncols=2;
ncols_write=10001;
all_h = [1, linspace(0,40,ncols_write-1)];
td = 2;

% S8 data
dbname = '/scratch/adrian/srvr_data_1.h5';
dsetname = '/lr15hr101.258639865h1T2/trials';
dsetname_decision = '/lr15hr101.258639865h1T2/decision_nonlin';
trial_data = h5read(dbname,dsetname,[1 1],[ncols ntrials]);
kappa = log(101.258639865/15);

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
