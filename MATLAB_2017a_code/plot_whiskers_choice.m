% plots whisker plots of choice matches
clear
rng('shuffle')
nsamples=100; % ideally 10K
for TN=100:100:200
    % nonlin
    % load modes
    load(['../data/mse_nonlin_fig4_iteration2_',num2str(TN),'trials.mat'])
    load(['../data/mse_lin_fig4_iteration2_',num2str(TN),'trials.mat'])

    for num_disc=1:500
        % set discounting parameter to mode value for the 4 fitted models
        
        % cross-validate on out-of-sample trials
        for num_sample=1:nsamples
            % generate clicks
            
            % run the 4 fitted models on this trial, and the 2 reference models
            
            
            % update the match counts variables
        end
    end
    lw=4;
    ms=8;
    fs=20;
    boxplot([accuracy_nonlinlin;...
        accuracy_nonlinnonlin;...
        accuracy_linlin;...
        accuracy_linnonlin]')
    set(findobj(gca,'type','line'),'linew',lw)
    set(gca,'linew',lw/2)
    ylim([.85,.9])
    ylabel('acc')
    ax=gca;ax.FontSize=fs;
    saveas(gcf, ['whisker_',num2str(TN),'.png'])
end

function acc = mapmode2match(modes, samples, sample_acc)
% map mode data to choice match proportion
%ARGS:
%   modes   - row vector of modes
%   sample_acc   - accuracy values corresponding to samples
%   samples - row vector of samples (discount. parameter values)
%RETURN: row vector of accuracy values
%NOTES: modes and samples must have same length
% samples must be sorted increasingly
% to get data do: load('../data/mse_nonlin_fig4_iteration2_100trials.mat')
if max(modes) > max(samples)
    error('max(modes)>max(samples)')
elseif min(modes) < min(samples)
    error('min(modes)<min(samples)')
else
    nmodes=length(modes);
    acc=zeros(size(modes));
    sorted_modes = sort(modes);
    for i=1:nmodes
        mode = sorted_modes(i);
        idx=find(mode<=samples,1);
        pidx=idx-1;
        if pidx > 0
            if abs(mode-samples(pidx)) < abs(mode-samples(idx))
                acc(i)=sample_acc(pidx);
            else
                acc(i)=sample_acc(idx);
            end
        end
    end
end
end