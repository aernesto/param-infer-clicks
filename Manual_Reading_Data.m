%% Description of data file
% the data files are located in '../synthetic_data/'
% and have an extension '.h5' which corresponds to the
% hdf5 data format. 

% let's first clear all our data
clear

% filename example: 
filename = '/home/adrian/programing/data/clicks/db1.h5';

%% Data structure 
% The data inside the .h5 file is stored inside a 'group'
% To find the group name, you can do the following
file_info = h5info(filename);
group_name = file_info.Groups.Name;

% The group itself contains 4 datasets: 'trials', 'trial_info',
% 'decision_lin' and 'decision_nonlin':
trials_dset_name=[group_name,'/trials'];
info_dset_name=[group_name,'/trial_info'];
nonlin_decision_dset_name=[group_name,'/decision_nonlin'];
lin_decision_dset_name=[group_name,'/decision_lin'];

% Uncomment the following to display information about the group on the 
% console:
% h5disp(filename, group_name)

%% Read off parameters
% The dataset named '/trials' contains click times and environment state
% values for all the trials in the group.
% To find out how many trials are in the file, we do:
trial_info = h5info(filename, trials_dset_name);
tot_num_trials = trial_info.Dataspace.Size(2);  % nb of trials in dataset

% To get the values of more parameters regarding the data, we access the 
% attributes of the other datasets, as follows:
h = h5readatt(filename, info_dset_name,'h');  % hazard rate
T = h5readatt(filename, info_dset_name,'T');  % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); 
high_rate = h5readatt(filename, info_dset_name,'high_click_rate'); 
S = h5readatt(filename, info_dset_name,'S');  % (high_rate-low_rate)/sqrt(low_rate+high_rate)
best_gamma = h5readatt(filename, lin_decision_dset_name, 'disc_lin');

%% Extract trials data
% To store all the click times and change point times in the group into a
% cell array, we do:
all_trials = h5read(filename, [group_name,'/trials']);
% row 1 in the cell corresponds to left trains
% row 2 in the cell corresponds to right trains
% row 3 in the cell corresponds to change point times
% There is one column in the cell per trial

% to get the start and end states of the environment for each trial, we do:
all_envt = h5read(filename, [group_name,'/trial_info'], [1 1], [2 Inf]);
% row 1 of the cell corresponds to starting state of the environment
% row 2 of the cell corresponds to ending state of the environment
% a value of 1 indicates H+, -1 indicates H-

% To get the left and right click trains of trial 50, for instance, as two
% column vectors, we do:
[left_train, right_train] = all_trials{1:2,50};  
% the correct choice on this trial H+ since:
all_envt(2,50) == 1

%% Extract decision data 
% The decision data of the linear (resp. nonlinear) model with the correct
% discounting parameter, for all trial, is given by the following row 
% vectors:
linear_decisions = h5read(filename,...
    [group_name,'/decision_lin'], [1 1], [1 Inf]);
nonlinear_decisions = h5read(filename,...
    [group_name,'/decision_nonlin'], [1 1], [1 Inf]);
% NOTE: 1 stands for H+ decision, -1 for H-, and 0 for an undecisive state.
% This happens when accumulation variable reaches end of trial at value 0
% (for instance if there is no click in the trial; or same number of left 
% and right clicks and discounting parameter is 0)


% As an example, accuracy of linear vs nonlinear model may be computed as 
% follows:
lin_acc = sum(linear_decisions == all_envt(2,:))/tot_num_trials
nonlin_acc = sum(nonlinear_decisions == all_envt(2,:))/tot_num_trials


