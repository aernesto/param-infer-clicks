%############ VALIDATED
% display('single epoch')
% tic
% v = end_point_nonlin(0.7, .5, 1);
% toc
%tic
%for i = 1:10000
%    v = end_point_nonlin(i/10000,i/9900,1);
%end
%toc

td = 2;
% ls = [1/3, pi/6, .8422343];
% rs = [1/6, 2/3, .75378];
kappa = log(101.258639865/15);
hh = 1;
% display('single trial')
% tic
% dd = decide1(td, ls, rs, kappa, hh, 0)
% toc

% hdf5 DB
dbname = '/scratch/adrian/sandbox_data.h5';
dsetname = '/lr15hr101.26h1T2/trials';
dsetname_decision = '/lr15hr101.26h1T2/decision_nonlin';
%h5disp(dbname, dsetname)

ntrials=90;
ncols=2;
ncols_write=101;
trial_data = h5read(dbname,dsetname,[1 1],[ncols ntrials]);
display('100 trials')

all_h = [1, linspace(0,40,ncols_write-1)];
tic
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
toc
