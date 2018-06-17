% test


% nonlinlin
% load modes
load('../data/mse_nonlin_fig4_iteration2_100trials.mat')
% load accuracy values
%load('../data/linacc.csv','-ascii')  % loads linacc
load('../data/linacc.csv','-ascii')  % loads nonlinacc
accuracy_nonlinlin = mapmode2acc(modes_nonlinlin, linspace(0,10,1000),...
    nonlinacc');
accuracy_nonlinnonlin = mapmode2acc(modes_nonlinnonlin,...
    linspace(0,10,1000), nonlinacc');

boxplot([accuracy_nonlinlin;accuracy_nonlinnonlin]')
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