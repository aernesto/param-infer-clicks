function yy = epn_AR(init_cond, end_time, hr)
    yy=zeros(size(init_cond));
    m=(init_cond ~= 0);
    yy(m) = 2*acoth(exp(2*hr(m)*end_time).*coth(init_cond(m)/2));
end
