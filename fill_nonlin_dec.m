function fill_nonlin_dec(lr, hr, kappa, h, td, dbname, ntrials, nsamples, sample_range)
% arguments passed to the script should be in the following order:
%    1. low rate
%    2. high rate
%    3. kappa
%    4. hazard rate
%    5. interrogation time
%    6. db filename
%    7. num_trials
%    8. num_samples
%    9. sample_range
%dbname = char(filename);
%dbname = filename;
%display(dbname)
ncols=2;
ncols_write=nsamples+1;
all_h = [h, linspace(sample_range(1),sample_range(end),ncols_write-1)];
ss=strip(strip(num2str(hr,'%.2f'),'0'),'.');
grp_name = ['/lr',num2str(lr),'hr',ss,'h',...
    num2str(h,'%.0f'), 'T',num2str(td, '%.0f')];
%display(grp_name)
dsetname = [grp_name,'/trials'];
dsetname_decision = [grp_name,'/decision_nonlin'];
trial_data = h5read(dbname,dsetname,[1 1],[ncols ntrials]);
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
h5writeatt(dbname,dsetname_decision,'true_h',h)
h5writeatt(dbname,dsetname_decision,'init_sample',sample_range(1))
h5writeatt(dbname,dsetname_decision,'end_sample',sample_range(2))
h5writeatt(dbname,dsetname_decision,'sample_step',diff(sample_range)/(ncols_write-2))
toc
