% checks convergence as fcn of nb of particles
clear
rng('shuffle')
lw=3; % line width for plots

hr = 25;
lr = 10;
T = .1; % trial duration
true_h = 1; % true hazard rate
k=log(hr/lr); % kappa for mean jump size in LLR at click

hs=3;%linspace(0,5,50)'; % values of h to try
ntrials=1;
lst = [.01]; 
rst = [.06]; 
parts = 500:100:2000;
noises = linspace(.5,1.5,3);
likelihoods = zeros(length(parts), length(noises));
tic
for kk = 1:length(noises)
    nsd=noises(kk);
    synthetic_decision = gauss_noise_nonlin_decide(T, lst, rst, k, true_h, 0, nsd);
    for jj = 1:length(parts)
        npart=parts(jj);            
        likelihoods(jj,kk)=lhd_nonlin_sing_tr_gauss_clicks(synthetic_decision, npart, lst, rst, T, k, hs, 0, nsd);
    end
    subplot(1,3,kk)
    plot(parts, likelihoods(:,kk),'LineWidth', 3)
    xlabel('nb of samples')
    ylabel('P(data|h)')
    title(['noise=',num2str(nsd)])
    ax=gca; ax.FontSize=16;
end
toc