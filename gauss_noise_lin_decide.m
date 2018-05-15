function decision = gauss_noise_lin_decide(left_train, right_train, gamma, kappa, noise_sd, seed, y0)
% returns decision, -1,0 or 1, using the linear model with Gaussian
% multiplicative noise applied to clicks
rng(seed);
noise_left_train = normrnd(kappa, noise_sd, size(left_train));
noise_right_train = normrnd(kappa, noise_sd, size(right_train));
yp=sum(noise_right_train.*exp(gamma*right_train));
ym=sum(noise_left_train.*exp(gamma*left_train));
decision=sign(y0+yp-ym);
end