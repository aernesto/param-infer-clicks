% analysis script
clear; clc;
%% specify database with choice data
dbname='~/programing/data/clicks/db1.h5';
% display parameters for current db
fprintf('disc and noise params below are for linear model\n')
p_lin=fetch_params(dbname,'lin')
fprintf('disc and noise params for nonlinear model are\n')
p_nonlin=fetch_params(dbname,'nonlin');
p_nonlin.disc
p_nonlin.noise

%% specify training and validation datasets
cutoff=.1; % percentage of trials to use for training
int_cutoff=floor(cutoff*p_lin.tot_db_trials);
fprintf('about to use %i trials for the fits\n',int_cutoff)
training_trials_range=[1,int_cutoff];
validation_trials_range=[int_cutoff+1,p_lin.tot_db_trials];

%% fit models (all four combinations)
types={'lin','nonlin'};

% specify parameters for the fits

% range of posterior support
lin_prior_range=[0,40];
nonlin_prior_range=[0,20];

% number of points to use for support of posterior
num_points_posterior_lin=600;
num_points_posterior_nonlin=600;

num_particles_lhd=400; % number of particles to use to estimate likelihood

point_estimates= containers.Map;

for i=1:2
    ref_type=types{i};
    for ii=1:2
        fitted_type=types{ii};
        [~,point_estimate]=fit_model(fitted_type,ref_type,...
            dbname,...
            lin_prior_range,...
            nonlin_prior_range,...
            training_trials_range,...
            num_points_posterior_lin,...
            num_points_posterior_nonlin,...
            num_particles_lhd);
        point_estimates([fitted_type,ref_type])=point_estimate;
    end
end

%% Compute analysis quantities
analysis_results=containers.Map;
for i=1:2
    ref_type=types{i};
    for ii=1:2
        fitted_type=types{ii};
        combination_string=[fitted_type,ref_type];
        fitted_params.disc=point_estimates(combination_string);
        
        % it is intentional that I set the noise from the fitted model
        % equal to the noise from the reference model.
        if strcmp(ref_type,'lin')
            fitted_params.noise=p_lin.noise;
            
            ref_params.disc=p_lin.disc;
            ref_params.noise=p_lin.noise;
            
        elseif strcmp(ref_type,'nonlin')
            fitted_params.noise=p_nonlin.noise;
            
            ref_params.disc=p_nonlin.disc;
            ref_params.noise=p_nonlin.noise;
        end
        
        % accuracy
        tmp_struct.acc=accuracy(fitted_type,fitted_params,dbname,...
            validation_trials_range);
        
        % predictive power 
        tmp_struct.pp=predictive_power(fitted_type, fitted_params,...
            ref_type, dbname, validation_trials_range);
        
        analysis_results(combination_string)=tmp_struct;
    end
    
    % get accuracy and predictive power of 'true model', which is the one that
    % was used to produce the validation data.
    
    % accuracy
    tmp_struct.acc=accuracy(ref_type,ref_params,dbname,...
        validation_trials_range);
    
    % predictive power
    tmp_struct.pp=predictive_power(ref_type, ref_params,...
        ref_type, dbname, validation_trials_range);
    
    analysis_results([ref_type,'true'])=tmp_struct;
end


%% display results
varNames={'L_acc','L_pp','NL_acc','NL_pp','T_acc','T_pp'};
rowNames={'ref_L','ref_NL'};
Results_Table=table('Size',[2,6],'VariableTypes',...
    {'double','double','double','double','double','double'},...
    'VariableNames',varNames,...
    'RowNames',rowNames);

% column for fitting linear model
Results_Table.L_acc=[analysis_results('linlin').acc;...
    analysis_results('linnonlin').acc];
Results_Table.L_pp=[analysis_results('linlin').pp;...
    analysis_results('linnonlin').pp];

% column for fitting nonlinear model
Results_Table.NL_acc=[analysis_results('nonlinlin').acc;...
    analysis_results('nonlinnonlin').acc];
Results_Table.NL_pp=[analysis_results('nonlinlin').pp;...
    analysis_results('nonlinnonlin').pp];

% column for true model
Results_Table.T_acc=[analysis_results('lintrue').acc;...
    analysis_results('linnonlin').acc];
Results_Table.T_pp=[analysis_results('nonlintrue').pp;...
    analysis_results('linnonlin').pp];

Results_Table

%% Auxiliary functions

% the following function might try to do too many things at the same time
function [posterior,estimate]=fit_model(fitted_type,ref_type,...
    dbname,...
    lin_prior_range,...
    nonlin_prior_range,...
    training_trials_range,...
    num_points_posterior_lin,...
    num_points_posterior_nonlin,...
    num_particles_lhd)
shuffle_db=false;
    if strcmp(fitted_type,'lin')
        [posterior,estimate]=fit_linear_model(ref_type, dbname,...
            lin_prior_range,training_trials_range,...
            num_points_posterior_lin,shuffle_db);
    elseif strcmp(fitted_type,'nonlin')
        [posterior, estimate]=fit_nonlinear_model(ref_type, dbname,...
            nonlin_prior_range,training_trials_range,...
            num_points_posterior_nonlin,num_particles_lhd,shuffle_db);
    end
end
