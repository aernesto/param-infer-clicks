function density = llh2density_AR(LLH, dh)
    vals = exp(LLH);
    N = dh * sum(vals);
    density = vals/N;
end
