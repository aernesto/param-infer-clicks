% plot whiskers of accuracy recovered from modes of posteriors for
% fits of stochastic models
clear
for TN=100:100:500  % loop over block size (in trial numbers)
    % nonlin
    % load modes
    load(['../data/mse_nonlin_fig4_iteration2_',num2str(TN),'trials.mat'])
    load(['../data/mse_lin_fig4_iteration2_',num2str(TN),'trials.mat'])
    % load accuracy values
    %load('../data/linacc.csv','-ascii')  % loads linacc
    load('../data/linacc_S3lr2.csv','-ascii')  % loads linacc for db S3lr2
    linacc=linacc_S3lr2;
    load('../data/nonlinacc.csv','-ascii')  % loads nonlinacc
    accuracy_nonlinlin = mapmode2acc(modes_nonlinlin, linspace(0,10,1000),...
        nonlinacc');
    accuracy_nonlinnonlin = mapmode2acc(modes_nonlinnonlin,...
        linspace(0,10,1000), nonlinacc');
    accuracy_linlin = mapmode2acc(modes_linlin, linspace(0,40,10000),...
        linacc');
    accuracy_linnonlin = mapmode2acc(modes_linnonlin,...
        linspace(0,40,10000), linacc');
    %g=['NLL','NLNL','LL','LNL'];%[ones(500,1),2*ones(500,1),3*ones(500,1),4*ones(500,1)];
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

function acc = mapmode2acc(modes, samples, sample_acc)
% map mode data to accuracy data
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