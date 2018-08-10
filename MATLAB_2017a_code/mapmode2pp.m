function PP = mapmode2pp(modes, samples, sample_pp)
% map mode data to PP data
%ARGS:
%   modes   - row vector of modes
%   sample_acc   - PP values corresponding to samples
%   samples - row vector of samples (discount. parameter values)
%RETURN: row vector of PP values
% samples must be sorted increasingly
% to get data do: load('../data/mse_nonlin_fig4_iteration2_100trials.mat')
if max(modes) > max(samples)
    warning('max(modes)>max(samples)')
elseif min(modes) < min(samples)
    error('min(modes)<min(samples)')
else
    nmodes=length(modes);
    PP=zeros(size(modes));
    sorted_modes = sort(modes);
    for i=1:nmodes
        mode = sorted_modes(i);
        idx=find(mode<=samples,1);
        pidx=idx-1;
        if pidx > 0
            if abs(mode-samples(pidx)) < abs(mode-samples(idx))
                PP(i)=sample_pp(pidx);
            else
                PP(i)=sample_pp(idx);
            end
        end
    end
end
end