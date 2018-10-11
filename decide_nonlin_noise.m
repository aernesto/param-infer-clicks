function decision = decide_nonlin_noise(trial_duration,...
    left_clicks, right_clicks, hh, init_cond, noise)
% runs N particles through the same trial with different noise
% realizations, where N is the number of distinct discounting parameter
% values to use. 
%
% So each particle has its own leak rate
%
% ARGS:
%   trial_duration: in sec
%   left_clicks:    col vec
%   right_clicks:   col vec
%   hh:             column vector of discounting parameter values to use
%   init_cond:      starting value of log-posterior odds ratio
%   noise:          2D noise matrix, dim1=click location in trial; 
%                   dim2=discounting value.
% RETURNS: a column vector of size size(hh) with choice values
% NOTES: rng('shuffle') should be called before function call
%   Called by: fill_nonlin_dec()
%   Calls: nonlin_eval() -- at end of script

    if size(noise,2)>1 && size(noise,1)==1
        noise=noise';
    end

    y = init_cond*ones(size(hh));
    t = 0;
    right_clicks_left = size(right_clicks, 1);
    left_clicks_left = size(left_clicks, 1);
    clicks_left = left_clicks_left + right_clicks_left;
    noise_idx=0;  % indexes click location
    while clicks_left > 0
        noise_idx = noise_idx+1;
        if right_clicks_left~=0 && left_clicks_left~=0
            if right_clicks(1) < left_clicks(1)
                nxt_click = right_clicks(1);
                right_clicks = right_clicks(2:end);
                right_clicks_left = right_clicks_left - 1;
                dwell = nxt_click - t;
                y = nonlin_eval(y, dwell, hh);
                y = y + noise(noise_idx,:)';
                t = nxt_click;
            elseif right_clicks(1) > left_clicks(1)
                nxt_click = left_clicks(1);
                left_clicks = left_clicks(2:end);
                left_clicks_left = left_clicks_left - 1;
                dwell = nxt_click - t;
                y = nonlin_eval(y, dwell, hh);
                y = y - noise(noise_idx,:)';
                t = nxt_click;
            else
                right_clicks_left = right_clicks_left - 1;
                left_clicks_left = left_clicks_left - 1;
                right_clicks = right_clicks(2:end);
                left_clicks = left_clicks(2:end);
            end
        elseif right_clicks_left~=0
            nxt_click = right_clicks(1);
            right_clicks = right_clicks(2:end);
            right_clicks_left = right_clicks_left - 1;
            dwell = nxt_click - t;
            y = nonlin_eval(y, dwell, hh);
            y = y + noise(noise_idx,:)';
            t = nxt_click;
        elseif left_clicks_left~=0
            nxt_click = left_clicks(1);
            left_clicks = left_clicks(2:end);
            left_clicks_left = left_clicks_left - 1;
            dwell = nxt_click - t;
            y = nonlin_eval(y, dwell, hh);
            y = y - noise(noise_idx,:)';
            t = nxt_click;
        end
        clicks_left = left_clicks_left + right_clicks_left;
    end
    dwell = trial_duration - t;
    y = nonlin_eval(y, dwell, hh);
    decision = sign(y);
end

function yy = nonlin_eval(init_cond, end_time, hr)
% evaluates evidence at the next click by solving nonlinear system.
%
% does this in vectorized form for several initial conditions and
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
