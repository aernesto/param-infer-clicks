function lklh=lhd_nonlin_sing_tr_gauss_clicks(dec_data,...
    npart, left_clicks, right_clicks, td, kappa, h_sample, init_cond,...
    noise_stdev, noise_bank)
% DESCR:
%
% ARGS:
%   dec_data: either -1 or 1, corresponds to value we want the lh of.
%   npart: number of particles
%   left_clicks:
%   right_clicks:
%   td: trial duration in sec
%   kappa:
%   h_sample: vector of values of h to use to compute decisions
%   init_cond:
%   noise_stdev: stdev of Gaussian multiplicative noise applied to clicks
%   noise_bank:
% RETURNS:

ncorrect = 0;
parfor part_nb = 1:npart
    dec = gauss_noise_nonlin_decide(td, left_clicks, right_clicks,...
        kappa, h_sample, init_cond, noise_stdev, noise_bank(part_nb,:)');
    if dec == dec_data
        ncorrect = ncorrect + 1;
    end
end
lklh = ncorrect / npart;
end