function [linlin, linnonlin, nonlinlin, nonlinnonlin]=count_perfect_samples(dbname, dsetname, numtrials, numsamples, init_g, last_g, init_h, last_h)
    % computes empirical probabilities of each sample being the correct one
    % returns a struct containing max probabilities and corresponding
    % sample values, for each model combination
    dsetname_decision_nonlin = [dsetname,'/decision_nonlin'];
    dsetname_decision_lin = [dsetname,'/decision_lin'];
    dsetname_info = [dsetname,'/trial_info'];
    true_param_lin = h5readatt(dbname,dsetname_decision_lin,'best_gamma');
    true_param_nonlin = h5readatt(dbname, dsetname_info, 'h');
    linear_decisions = h5read(dbname, dsetname_decision_lin, [2,1], [numsamples, numtrials]);
    nonlinear_decisions = h5read(dbname, dsetname_decision_nonlin, [2,1], [numsamples, numtrials]);
    reference_decision_linear = h5read(dbname, dsetname_decision_lin, [1,1], [1, numtrials]);
    reference_decision_nonlinear = h5read(dbname, dsetname_decision_nonlin, [1,1], [1, numtrials]);
    lin_param_values = linspace(init_g,last_g,numsamples);
    nonlin_param_values = linspace(init_h, last_h, numsamples);
    % linlin
    bool_lin_lin = (linear_decisions == reference_decision_linear);
    proba_linear = mean(bool_lin_lin, 2);
    linlin.maxprob = max(proba_linear);
    linlin.count = sum(proba_linear==linlin.maxprob);
    linlin.values = lin_param_values(proba_linear==linlin.maxprob);
    linlin.true_param_dec = true_param_lin;
    linlin.true_param_fit = true_param_lin;
    % nonlinnonlin
    bool_nonlin_nonlin = (nonlinear_decisions == reference_decision_nonlinear);
    proba_nonlinear = mean(bool_nonlin_nonlin, 2);
    nonlinnonlin.maxprob = max(proba_nonlinear);
    nonlinnonlin.count = sum(proba_nonlinear==nonlinnonlin.maxprob);
    nonlinnonlin.values = nonlin_param_values(proba_nonlinear==nonlinnonlin.maxprob);
    nonlinnonlin.true_param_dec = true_param_nonlin;
    nonlinnonlin.true_param_fit = true_param_nonlin;
    % linnonlin - fit linear model to nonlinear data
    bool_lin_nonlin = (linear_decisions == reference_decision_nonlinear);
    proba_linnonlinear = mean(bool_lin_nonlin, 2);
    linnonlin.maxprob = max(proba_linnonlinear);
    linnonlin.count = sum(proba_linnonlinear==linnonlin.maxprob);
    linnonlin.values = lin_param_values(proba_linnonlinear==linnonlin.maxprob);
    linnonlin.true_param_dec = true_param_nonlin;
    linnonlin.true_param_fit = true_param_lin;
    % nonlinlin - fit nonlinear model to linear data
    bool_nonlin_lin = (nonlinear_decisions == reference_decision_linear);
    proba_nonlinlinear = mean(bool_nonlin_lin, 2);
    nonlinlin.maxprob = max(proba_nonlinlinear);
    nonlinlin.count = sum(proba_nonlinlinear==nonlinlin.maxprob);
    nonlinlin.values = nonlin_param_values(proba_nonlinlinear==nonlinlin.maxprob);
    nonlinlin.true_param_dec=true_param_lin;
    nonlinlin.true_param_fit=true_param_nonlin;
end
