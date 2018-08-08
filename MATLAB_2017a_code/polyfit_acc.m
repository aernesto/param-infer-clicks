% plot acc of L and NL for range of discounting parameters
clear
NL_data=load('../data/acc_NL_range_h.mat');
L_data=load('../data/acc_L_corrected_1.mat');
%a2=load('../data/acc_L_corrected_2.mat'); % alternative data

% several plotting possibilities
xNL=NL_data.hs; yNL=NL_data.acc;
xL=L_data.gammas; yL=L_data.acc;

yyaxis left
plot(yNL,xNL)
yyaxis right
plot(yL,xL)

% subplot(2,1,1)
% plot(NL_data.hs,NL_data.acc)
% xlabel('h')
% subplot(2,1,2)
% plot(L_data.gammas,L_data.acc)
% subplot(3,1,3)

% plot(a2.gammas,a2.acc)
% plot(a.gammas,a.acc,a2.gammas,a2.acc)
% plot(a.gammas,smooth(a.acc),a2.gammas,smooth(a2.acc))

p = polyfit(L_data.gammas,L_data.acc,7);

% plot data and fitted polynomial on top
plot(L_data.gammas,polyval(p,L_data.gammas),L_data.gammas,L_data.acc)

% max acc of noisy data
max(L_data.acc)
% ans =
% 
%     0.8344

% max acc of fitted polyomial
max(polyval(p,L_data.gammas))
% ans =
% 
%     0.8340

[m1,m2]=max(polyval(p,L_data.gammas))

% m1 =
% 
%     0.8340
% 
% 
% m2 =
% 
%    114

[n1,n2]=max(L_data.acc)

% n1 =
% 
%     0.8344
% n2 =
% 
%    114

% argmax for L is 114 
L_data.gammas(114)

% ans =
% 
%     5.6500

% p2 = polyfit(a2.gammas,a2.acc,7);
% plot(a.gammas,polyval(p,a.gammas),a2.gammas,polyval(p2,a2.gammas))
% [o1,o2]=max(polyval(p2,a2.gammas))
% 
% o1 =
% 
%     0.8340
% 
% 
% o2 =
% 
%    114

%--------------------- NL-------------------------

p3 = polyfit(NL_data.hs,NL_data.acc,10);
plot(NL_data.hs,NL_data.acc,NL_data.hs,polyval(p3,NL_data.hs))
[b1,b2]=max(NL_data.acc);
plot(NL_data.hs,NL_data.acc,[NL_data.hs(b2),NL_data.hs(b2)],[0.6,b1])
NL_data.hs(b2)

% ans =
% 
%     0.4500
plot(NL_data.hs/NL_data.hs(46),NL_data.acc,L_data.gammas/L_data.gammas(114),L_data.acc)
legend('NL','L')
