function heights_mat=hist_part_filt(npart, left_clicks, right_clicks, td,...
    g_sample, init_cond, noise_bank)
% Compute and return the vectors of landing heights of particles with the 
% linear model.
% ARGS:
%   npart:      number of independent particles to use for each param value
%   left_clicks:
%   right_clicks:
%   td:         trial duration in sec
%   kappa:      jump size in evidence at click times
%   h_sample:   vector of values of h to use to compute decisions
%   init_cond:  
%   noise_stdev:stdev of Gaussian multiplicative noise applied to clicks
%   noise_bank: 3D matrix of i.i.d Gaussians to use for click noise
%
% RETURNS: matrix with col vectors = landing heights of each particle
% NOTES: 
%   Called by: histogram_test.m

heights_mat = zeros(npart,length(g_sample));
for part_nb = 1:npart
    heights_mat(part_nb,:) = ev_height_lin(td, left_clicks,...
        right_clicks, g_sample',...
        init_cond, squeeze(noise_bank(part_nb,:,:)));
end
end
