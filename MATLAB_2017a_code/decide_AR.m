function decision = decide_AR(trial_duration,...
    left_clicks, right_clicks, ~, hh, init_cond, ~, noise)
% runs N particles through the same trial with different noise
% realizations, where N is the number of distinct discounting parameter
% values to use. 
%
% So each particle has its own leak rate
%
% ARGS:
%   trial_duration: 
%   left_clicks:
%   right_clicks:
%   hh:             column vector of discounting parameter values to use
%   init_cond:
%   noise:          2D noise matrix, dim1=click location in trial; 
%                   dim2=discounting value.
% RETURNS: a column vector of size size(hh)
% NOTES: rng('shuffle') should be called before function call
%   Called by: lhd_AR.m; lin_stoch_fit.m

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
                y = epn_AR(y, dwell, hh);
                y = y + noise(noise_idx,:)';
                t = nxt_click;
            elseif right_clicks(1) > left_clicks(1)
                nxt_click = left_clicks(1);
                left_clicks = left_clicks(2:end);
                left_clicks_left = left_clicks_left - 1;
                dwell = nxt_click - t;
                y = epn_AR(y, dwell, hh);
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
            y = epn_AR(y, dwell, hh);
            y = y + noise(noise_idx,:)';
            t = nxt_click;
        elseif left_clicks_left~=0
            nxt_click = left_clicks(1);
            left_clicks = left_clicks(2:end);
            left_clicks_left = left_clicks_left - 1;
            dwell = nxt_click - t;
            y = epn_AR(y, dwell, hh);
            y = y - noise(noise_idx,:)';
            t = nxt_click;
        end
        clicks_left = left_clicks_left + right_clicks_left;
    end
    dwell = trial_duration - t;
    y = epn_AR(y, dwell, hh);
    decision = sign(y);
end
