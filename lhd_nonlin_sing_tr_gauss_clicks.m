function lklh=lhd_nonlin_sing_tr_gauss_clicks(dec_data, npart,...
    left_clicks, right_clicks, td, h_sample, init_cond, noise_bank)
% Estimates a vector of likelihood values, one per discounting parameter 
% value. Each likelihood is estimated via sampling (particle filter).
%
% The likelihood of a value is the probability that the model makes the 
% same decision as the data when using this discounting value for this 
% trial. 
%
% ARGS:
%   dec_data:   either -1 or 1, corresponds to value we want the lh of.
%   npart:      number of independent particles to use for each param value
%   left_clicks:
%   right_clicks:
%   td:         trial duration in sec
%   h_sample:   vector of values of h to use to compute decisions
%   init_cond:  
%   noise_bank: 3D matrix of i.i.d Gaussians to use for click noise
%
% RETURNS: col vector of likelihoods
% NOTES: 
%   Called by: fit_nonlinear_model()
%   Calls: decide_nonlin_noise()

ncorrect = zeros(size(h_sample));
for part_nb = 1:npart
    dec = decide_nonlin_noise(td, left_clicks, right_clicks, h_sample,...
        init_cond, squeeze(noise_bank(part_nb,:,:)));
    ncorrect(dec == dec_data) = ncorrect(dec == dec_data) + 1;
end
lklh = ncorrect / npart;
end