function heights = ev_height_lin(trial_duration,...
    left_clicks, right_clicks, hh, init_cond, noise)
% runs N particles through the same trial with different noise
% realizations, where N is the number of distinct discounting parameter
% values to use. 
%
% So each particle has its own leak rate
%
% ARGS:
%   trial_duration: scalar
%   left_clicks:    col vec of click times
%   right_clicks:   col vec of click times
%   hh:             column vector of discounting parameter values to use
%   init_cond:      scalar
%   noise:          2D noise matrix, dim1=click location in trial; 
%                   dim2=discounting value.
% RETURNS: a row vector of landing heights of the particles
% NOTES:
%   THIS FUNCTION NEEDS TESTING!
%   rng('shuffle') should be called before function call
%   Called by: hist_part_filt.m

% size(trial_duration)
% size(left_clicks)
% size(right_clicks)
% size(hh)
% size(init_cond)
% size(noise)
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
                y = linev(y, dwell, hh);
                y = y + noise(noise_idx,:)';
                t = nxt_click;
            elseif right_clicks(1) > left_clicks(1)
                nxt_click = left_clicks(1);
                left_clicks = left_clicks(2:end);
                left_clicks_left = left_clicks_left - 1;
                dwell = nxt_click - t;
                y = linev(y, dwell, hh);
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
            y = linev(y, dwell, hh);
            y = y + noise(noise_idx,:)';
            t = nxt_click;
        elseif left_clicks_left~=0
            nxt_click = left_clicks(1);
            left_clicks = left_clicks(2:end);
            left_clicks_left = left_clicks_left - 1;
            dwell = nxt_click - t;
            y = linev(y, dwell, hh);
            y = y - noise(noise_idx,:)';
            t = nxt_click;
        end
        clicks_left = left_clicks_left + right_clicks_left;
    end
    dwell = trial_duration - t;
    y = linev(y, dwell, hh);
    heights = y';
end
function k=linev(init_cond, end_time, hr)
% evaluates evidence at the next click by solving linear system.
%
% does this in vectorized form form several initial conditions and
% discounting parameter values
%
% ARGS:
%   init_cond: column vector of initial values for evidence variable
%   end_time:  scalar
%   hr:        column vector of discounting parameter values to use
% RETURNS: column vector of evidence values
    k = init_cond.*exp(-hr*end_time);
end
