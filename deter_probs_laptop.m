% computes empirical probabilities of each sample being the correct one
clear
% explore data
dbname = 'data/small_data_TEST.h5';
%dsetname = '/lr1hr11.68h1T2/trials';
dsetname_decision_nonlin = '/lr1hr11.68h1T2/decision_nonlin';
dsetname_decision_lin = '/lr1hr11.68h1T2/decision_lin';
%h5disp(dbname)
numtrials=1000;
numsamples = 100;
linear_decisions = h5read(dbname, dsetname_decision_lin, [2,1], [numsamples, numtrials]);
nonlinear_decisions = h5read(dbname, dsetname_decision_nonlin, [2,1], [numsamples, numtrials]);

%size(linear_decisions)
%size(nonlinear_decisions)
reference_decision_linear = h5read(dbname, dsetname_decision_lin, [1,1], [1, numtrials]);
reference_decision_nonlinear = h5read(dbname, dsetname_decision_nonlin, [1,1], [1, numtrials]);

true_param_lin = h5readatt(dbname,dsetname_decision_lin,'best_gamma');
true_param_nonlin = 1;
param_values = linspace(0,40,numsamples);

% linlin
bool_lin_lin = (linear_decisions == reference_decision_linear);
proba_linear = mean(bool_lin_lin, 2);
figure(1)
plot(param_values,proba_linear, 'LineWidth', 3)
hold on
plot([true_param_lin true_param_lin],[0 1],'red','LineWidth', 2)
hold off
title('Lin fit to lin')
xlabel('Disc. param. gamma')
ylabel('% compatible')
legend('proportion trials','true gamma (best)')
ax=gca;
ax.FontSize=20;

% nonlinnonlin
bool_nonlin_nonlin = (nonlinear_decisions == reference_decision_nonlinear);
proba_nonlinear = mean(bool_nonlin_nonlin, 2);
figure(2)
plot(param_values,proba_nonlinear, 'LineWidth', 3)
hold on
plot([true_param_nonlin true_param_nonlin],[0 1],'red','LineWidth', 2)
hold off
title('Nonlin fit to nonlin')
xlabel('Discount. param. h')
ylabel('% compatible')
legend('proportion trials','true hazard rate')
ax=gca;
ax.FontSize=20;

% linnonlin - fit linear model to nonlinear data
bool_lin_nonlin = (linear_decisions == reference_decision_nonlinear);
proba_linnonlinear = mean(bool_lin_nonlin, 2);

figure(3)
plot(param_values,proba_linnonlinear, 'LineWidth', 3)
hold on
plot([true_param_lin true_param_lin],[0 1],'red','LineWidth', 2)
hold off
title('Lin fit to nonlin')
xlabel('Disc. param. gamma')
ylabel('% compatible')
legend('proportion trials','best gamma')
ax=gca;
ax.FontSize=20;

% nonlinlin - fit nonlinear model to linear data
bool_nonlin_lin = (nonlinear_decisions == reference_decision_linear);
proba_nonlinlinear = mean(bool_nonlin_lin, 2);
figure(4)
plot(param_values,proba_nonlinlinear,'LineWidth', 3)
hold on
plot([true_param_nonlin true_param_nonlin],[0 1],'red','LineWidth', 2)
hold off
title('Nonlin fit to lin')
xlabel('Discount. param. h')
ylabel('% compatible')
legend('proportion trials','true h')
ax=gca;
ax.FontSize=20;
