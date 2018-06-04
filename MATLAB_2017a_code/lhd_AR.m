function lklh=lhd_AR(dec_data, npart, left_clicks, right_clicks, td, kappa, h_sample, init_cond, noise_stdev, noise_bank)
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
ncorrect = zeros(size(h_sample));
for part_nb = 1:npart
    dec = decide_AR(td, left_clicks, right_clicks, kappa, h_sample, init_cond, noise_stdev, noise_bank(part_nb,:));
    ncorrect(dec == dec_data) = ncorrect(dec == dec_data) + 1;
end
lklh = ncorrect / npart;
end
