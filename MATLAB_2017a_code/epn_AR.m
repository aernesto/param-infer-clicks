function yy = epn_AR(init_cond, end_time, hr)
% evaluates evidence at the next click by solving nonlinear system.
%
% does this in vectorized form form several initial conditions and
% discounting parameter values
%
% ARGS:
%   init_cond: column vector of initial values for evidence variable
%   end_time:  scalar
%   hr:        column vector of discounting parameter values to use
% RETURNS: column vector of evidence values
% NOTES: 
%   Called by decide_AR.m
    yy=zeros(size(init_cond));
    % we can only use the closed-form solution when evidence is nonzero
    m=(init_cond ~= 0);
    try
        yy(m) = 2*acoth(exp(2*hr(m)*end_time).*coth(init_cond(m)/2));
    catch ME
        size(yy)
        size(m)
        size(hr)
        size(init_cond)
        rethrow(ME)
    end
end
