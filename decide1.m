function decision = decide1(trial_duration, left_clicks, right_clicks, kappa, hh, init_cond)
    y = init_cond;
    t = 0;
    right_clicks_left = size(right_clicks, 1);
    left_clicks_left = size(left_clicks, 1);
    clicks_left = left_clicks_left + right_clicks_left;
    while clicks_left > 0
        if right_clicks_left~=0 && left_clicks_left~=0
            if right_clicks(1) < left_clicks(1)
                nxt_click = right_clicks(1);
                right_clicks = right_clicks(2:end);
                right_clicks_left = right_clicks_left - 1;
                dwell = nxt_click - t;
                y = end_point_nonlin(y, dwell, hh);
                y = y + kappa;
                t = nxt_click;
            elseif right_clicks(1) > left_clicks(1)
                nxt_click = left_clicks(1);
                left_clicks = left_clicks(2:end);
                left_clicks_left = left_clicks_left - 1;
                dwell = nxt_click - t;
                y = end_point_nonlin(y, dwell, hh);
                y = y - kappa;
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
            y = end_point_nonlin(y, dwell, hh);
            y = y + kappa;
            t = nxt_click;
        elseif left_clicks_left~=0
            nxt_click = left_clicks(1);
            left_clicks = left_clicks(2:end);
            left_clicks_left = left_clicks_left - 1;
            dwell = nxt_click - t;
            y = end_point_nonlin(y, dwell, hh);
            y = y - kappa;
            t = nxt_click;
        end
        clicks_left = left_clicks_left + right_clicks_left;
    end
    dwell = trial_duration - t;
    y = end_point_nonlin(y, dwell, hh);
    decision = sign(y);
end
