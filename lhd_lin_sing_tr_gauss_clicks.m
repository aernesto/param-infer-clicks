function lklh = lhd_lin_sing_tr_gauss_clicks(dec, noise_stdev, kappa, T, left_clicks, right_clicks, gammas)
% ARGS:
% dec = 1 or -1
% clicks must be row vectors (one of them at most may be empty)
% gammas must be a col vector or scalar
% RETURNS:
% returns NaN if both click trains are empty, otherwise, returns column
% vector with log probabilities, one entry per gamma value.
% DESCR:
% computes the log-likelihood of choosing dec (+1 or -1), given the click trains,
% the discounting rate gamma, and the stdev of the Gaussian noise applied
% to each click
if size(gammas,2) > 1
    error('gammas should be a column vector')
elseif isempty([right_clicks,left_clicks]) 
    error('both click trains are empty')
elseif isempty(right_clicks) 
    pos1 = sum(exp(gammas*reshape(right_clicks,1,[])),2);
    pos2 = sum(exp(2*gammas*reshape(right_clicks,1,[])),2);
    neg1 = sum(exp(gammas*left_clicks),2);
    neg2 = sum(exp(2*gammas*left_clicks),2);
elseif isempty(left_clicks) 
    pos1 = sum(exp(gammas*right_clicks),2);
    pos2 = sum(exp(2*gammas*right_clicks),2);
    neg1 = sum(exp(gammas*reshape(left_clicks,1,[])),2);
    neg2 = sum(exp(2*gammas*reshape(left_clicks,1,[])),2);
else
    pos1 = sum(exp(gammas*right_clicks),2);
    pos2 = sum(exp(2*gammas*right_clicks),2);
    neg1 = sum(exp(gammas*left_clicks),2);
    neg2 = sum(exp(2*gammas*left_clicks),2);
end
meanTerm = single(dec) * kappa * exp(-gammas*T) .* (pos1 - neg1);
varTerm = noise_stdev^2 * exp(-2*gammas*T) .* (pos2 + neg2);
lklh = log(normcdf(meanTerm ./ sqrt(varTerm)));
end