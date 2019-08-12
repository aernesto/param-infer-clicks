clear
ntrials=1000;
rng('shuffle')
ndiscount=200; % number of discounting parameter values to try
hstart=0;hend=10; % range should be large enough for normalization
hs=linspace(hstart,hend,ndiscount)'; % values of h to try
dh=(hend-hstart)/(ndiscount-1);
filename = '/home/adrian/Git/GitHub/work/param-infer-clicks/data/S3lr5h1T2tr10000sp1000.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
info_dset_name=[group_name,'/trial_info'];
true_h = h5readatt(filename, info_dset_name,'h');  % hazard rate
T = h5readatt(filename, info_dset_name,'T');  % Trial duration in sec
low_rate = h5readatt(filename, info_dset_name,'low_click_rate'); 
high_rate = h5readatt(filename, info_dset_name,'high_click_rate');
true_gamma = h5readatt(filename, [group_name,'/decision_lin'],'best_gamma');
k=log(high_rate/low_rate);
all_trials = h5read(filename, [group_name,'/trials']);
tot_trials_db = size(all_trials,2);
all_trials = all_trials(1:2,:);
npart = 800;
nsd=1; % Gaussian noise applied to click height

nruns=100;
mse_nonlinnonlin=0; mse_nonlinlin=0;
% store the modes of posteriors
modes_nonlinlin=zeros(1,nruns); modes_nonlinnonlin=modes_nonlinlin;
tic
for run=1:nruns
    % shuffle trial order
    all_trials = all_trials(:,randperm(tot_trials_db));
    llh = zeros(ndiscount,2); %col1=nonlin synthetic dec;col2=lin synth dec
    parfor trn=1:ntrials
        [lst, rst]=all_trials{:,trn};
        total_clicks = length(lst)+length(rst);
    
        % synthetic decision computed with linear model
        synthetic_decision_lin = gauss_noise_lin_decide(lst, rst,...
            true_gamma, k, nsd, 0);
        % compute reference decision with nonlinear model
        synthetic_decision_nonlin = decide_AR(T,...
            lst, rst, NaN, true_h, 0, NaN, normrnd(k, nsd, [total_clicks, 1]));
        % flip a coin if either decision was 0
        if synthetic_decision_lin == 0
            synthetic_decision_lin = sign(rand-0.5);
        end
        if synthetic_decision_nonlin == 0
            synthetic_decision_nonlin = sign(rand-0.5);
        end
    
        % generate upfront all the Gaussian r.v. needed
        Gaussian_bank_nonlinlin = normrnd(k, nsd, [npart, total_clicks, ndiscount]);
        Gaussian_bank_nonlinnonlin = normrnd(k, nsd, [npart, total_clicks, ndiscount]);
        
        % likelihood for nonlinnonlin
        lhv_nonlinnonlin=lhd_AR(synthetic_decision_nonlin, npart, lst, rst, T, k,...
            hs, 0, nsd, Gaussian_bank_nonlinnonlin);
        lhv_nonlinnonlin(lhv_nonlinnonlin<eps) = eps;
        
        % likelihood for nonlinlin
        lhv_nonlinlin=lhd_AR(synthetic_decision_lin, npart, lst, rst, T, k,...
            hs, 0, nsd, Gaussian_bank_nonlinlin);
        % prevent 0 likelihood
        lhv_nonlinlin(lhv_nonlinlin<eps) = eps;
        
        llh=llh+[log(lhv_nonlinnonlin), log(lhv_nonlinlin)];
    end
    
    % shift log-likelihood up to avoid numerical errors
    [max_nonlinlin,idx1]=max(llh(:,2));
    modes_nonlinlin(run)=hs(idx1);
    [max_nonlinnonlin,idx2]=max(llh(:,1));
    modes_nonlinnonlin(run)=hs(idx2);
    llh(:,1)=llh(:,1)+abs(max_nonlinnonlin);
    llh(:,2)=llh(:,2)+abs(max_nonlinlin);
    
    density_nonlinnonlin=llh2density_AR(llh(:,1),dh);
    mse_nonlinnonlin=mse_nonlinnonlin+dh*sum(((hs-true_h).^2).*density_nonlinnonlin);
    density_nonlinlin=llh2density_AR(llh(:,2),dh);
    mse_nonlinlin=mse_nonlinlin+dh*sum(((hs-true_h).^2).*density_nonlinlin);
end
mse_nonlinnonlin=mse_nonlinnonlin/nruns;
mse_nonlinlin=mse_nonlinlin/nruns;
toc
fname=['mse_nonlin_fig4_iteration3_',num2str(ntrials),'trials'];
save(['/home/adrian/programing/data/clicks/fits/',fname,'.mat'],'mse_nonlinnonlin',...
      'mse_nonlinlin', 'modes_nonlinlin', 'modes_nonlinnonlin')
