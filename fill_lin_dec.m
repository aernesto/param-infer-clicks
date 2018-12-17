function fill_lin_dec(disc, dbname, sigma)
% arguments passed to the script should be in the following order:
%    1. decision_hazard=hazard rate used by observer model to discount evidence
%    2. dbname=db filename (full path)
%    3. sigma=noise width on clicks height (stdev)
% NOTE: This function calls 'decide_nonlin_noise()' function
%       Also, there might be a way to parallelize the code, if speed is an
%       issue

% number of rows to read off for the stimulus data
ncols=2; 

% number of rows to write to db
ncols_write=1;

% fetch data from database
file_info = h5info(dbname);
grp_name = file_info.Groups.Name;
dsetname = [grp_name,'/trials'];
dsetname_decision = [grp_name,'/decision_lin'];
info_dset_name=[grp_name,'/trial_info'];

trial_info = h5info(dbname, dsetname);
ntrials = trial_info.Dataspace.Size(2);  % nb of trials in dataset

td = h5readatt(dbname, info_dset_name,'T');  % Trial duration in sec
lr = h5readatt(dbname, info_dset_name,'low_click_rate'); 
hr = h5readatt(dbname, info_dset_name,'high_click_rate'); 

    % stimulus data
trial_data = h5read(dbname,dsetname,[1 1],[ncols ntrials]);

% clicks mean height
kappa=log(hr/lr);

tic

% turn off warning generated by h5write call below
warning('off','MATLAB:imagesci:hdf5dataset:datatypeOutOfRange')

for i = 1:ntrials
    ls = trial_data{1,i}; % col vector
    rs = trial_data{2,i};

    dec = zeros(ncols_write,1);
    parfor tn=1:ncols_write
        % careful with rng here.
        rng('shuffle');
        dec(tn) = decide_lin_noise(ls, rs, disc, kappa, sigma, 0);
    end
    start=[1 i];
    count=[ncols_write 1];
    
    % write single decision datum to db
    h5write(dbname, dsetname_decision, dec, start, count)
end

h5writeatt(dbname,dsetname_decision,'disc_lin',disc)
h5writeatt(dbname, dsetname_decision,'noise',sigma)
toc
end