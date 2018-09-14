function density = llh2density_AR(LLH, dh)
% DESCR: converts vector of log likelihood values to a vector of density
% values
% ARGS: 
%   LLH - col vector of log likelihoods 
%   dh  - scalar
% RETURNS: col vector of pdf values, same length as LLH
% NOTES:
%   currently I am not sure how I should deal with very negative LLH values
    vals = exp(LLH);
    if max(vals) < eps
        density = zeros(size(LLH));
    else
        N = dh * sum(vals);
        density = vals/N;
    end
end
