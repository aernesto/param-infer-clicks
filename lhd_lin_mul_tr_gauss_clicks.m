function posterior = lhd_lin_mul_tr_gauss_clicks(dbname, grp_name, ntrials, noise_stdev, gammas)
ncols=2;
dsetname = [grp_name,'/trials'];
info_dset = [grp_name,'/trial_info'];
dec_dset = [grp_name,'/decision_lin'];

hr = h5readatt(dbname, info_dset, 'high_click_rate');
lr = h5readatt(dbname, info_dset, 'low_click_rate');
T = h5readatt(dbname, info_dset, 'T');

kappa = log(hr/lr);

trial_data = h5read(dbname, dsetname, [1 1], [ncols ntrials]);
dec_data = h5read(dbname, dec_dset, [1 1], [1 ntrials]);
%dec_data(1)
posterior = ones(size(gammas))/(gammas(end)-gammas(1));
tic
for i = 1:ntrials
    lst = trial_data{1,i}; % col vector
    rst = trial_data{2,i}; % col vector
    try
        posterior = posterior .* lhd_lin_sing_tr_gauss_clicks(dec_data(i), noise_stdev,...
            kappa, T, lst', rst', gammas);
    catch ME
        rethrow(ME)
        %warning(['two empty trains in trial ', num2str(i)]);
    end
end
toc
end