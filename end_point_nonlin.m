function yy = end_point_nonlin(init_cond, end_time, hr)
%DESCR: evaluates evidence at the next click by solving nonlinear system
%
%ARGS:
%   init_cond: 
%   end_time:
%   hr: 
%RETURNS:
    if init_cond == 0
        yy = 0;
    else
        yy = 2*acoth(exp(2*hr*end_time)*coth(init_cond/2));
    end
end
