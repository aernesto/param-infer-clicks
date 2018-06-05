% bug example for linear model fit:
clear
% if click times turned into single type with single(), then bug appears
% create right train
right_train=[0.2073,0.8918,1.1884,1.5134,1.6278,1.6789,1.7848,1.9237,1.9770];
% create left train
left_train=[0.2237,0.4118,0.4920,0.7067,0.7512,0.8246,0.8424,0.9041,...
            1.0474,1.2204,1.2526, 1.2910,1.4266,1.5040,1.5906,1.9554];

kappa=1.38629436111989057245;
decision_data=1;
noise_stdev= 1.5;
trial_duration=2;
discounting_rates=linspace(0,40,800)';
log_likelihoods = lhd_lin_sing_tr_gauss_clicks(decision_data, noise_stdev, kappa,...
    trial_duration, left_train, right_train, discounting_rates);
plot(log_likelihoods)