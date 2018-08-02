% compute predictive power NL-NL according to methodology described below:
% https://paper.dropbox.com/doc/Stochastic-model-fitting-URuXW5PKABizdbMJB4Bkb#:uid=862179754392750879685852&h2=Definition-of-predictive-power

% Only focus on NL-NL, with m_f and m_d using true discounting param

clear

%1. -----------bin q------------------------------------------------------%
%load('../data/nonlin_q_density_10000trials.mat')
load('../data/nonlin_q_density_noise_4_10000trials.mat')
nbins=hist.NumBins;  % actual number of bins for histogram and density
qbins=hist.BinEdges;  % really bin edges as needed by histogram()
qvalues=zeros(1,nbins); % mid-points of bins
for i=1:nbins
    qvalues(i)=(qbins(i+1)+qbins(i))/2;
end
    
%2. -----------get density of q across trials-----------------------------%

q_density=hist.Values;

%3. --------get probability of agreement as fcn of q----------------------%

agreement=qvalues.^2+(1-qvalues).^2;

%4. ----------------compute predictive power------------------------------%

bin_width=1/nbins; % equivalently bin_width=hist.BinWidth;

PP=bin_width*sum(agreement.*q_density)