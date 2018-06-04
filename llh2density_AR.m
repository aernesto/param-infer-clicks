function density = llh2density_AR(LLH, dh)
% DESCR: converts vector of log likelihood values to a vector of density
% values
% ARGS: LLH ; dh
% RETURNS:
% NOTES:
    vals = exp(LLH);
    N = dh * sum(vals);
    density = vals/N;
end
