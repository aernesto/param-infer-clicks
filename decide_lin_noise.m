function decision = decide_lin_noise(left_train, right_train,...
    gamma, kappa, noise_sd, y0)
% returns decision, -1,0 or 1, using the linear model with Gaussian
% multiplicative noise applied to clicks
% NOTES:
%   Called by: fill_lin_dec()
%   rng('shuffle') must be called before the function
    noise_left_train = normrnd(kappa, noise_sd, size(left_train));
    noise_right_train = normrnd(kappa, noise_sd, size(right_train));
    yp=sum(noise_right_train.*exp(gamma*right_train));
    ym=sum(noise_left_train.*exp(gamma*left_train));
    decision=sign(y0+yp-ym);
end
