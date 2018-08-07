b=load('../data/acc_NL_range_h.mat');
a=load('../data/acc_L_corrected_1.mat');
%a2=load('../data/acc_L_corrected_2.mat'); % alternative data
subplot(3,1,1)
plot(b.hs,b.acc)
xlabel('h')
subplot(3,1,2)
plot(a.gammas,a.acc)
subplot(3,1,3)
% plot(a2.gammas,a2.acc)
% plot(a.gammas,a.acc,a2.gammas,a2.acc)
% plot(a.gammas,smooth(a.acc),a2.gammas,smooth(a2.acc))
p = polyfit(a.gammas,a.acc,7);
plot(a.gammas,polyval(p,a.gammas),a.gammas,a.acc)
max(a.acc)

% ans =
% 
%     0.8344

max(polyval(p,a.gammas))

% ans =
% 
%     0.8340

[m1,m2]=max(polyval(p,a.gammas))

% m1 =
% 
%     0.8340
% 
% 
% m2 =
% 
%    114

[n1,n2]=max(a.acc)

% n1 =
% 
%     0.8344
% n2 =
% 
%    114
% 
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

p3 = polyfit(b.hs,b.acc,10);
plot(b.hs,b.acc,b.hs,polyval(p3,b.hs))
[b1,b2]=max(b.acc);
plot(b.hs,b.acc,[b.hs(b2),b.hs(b2)],[0.6,b1])
b.hs(b2)

% ans =
% 
%     0.4500
plot(b.hs/b.hs(46),b.acc,a.gammas/a.gammas(114),a.acc)
legend('NL','L')
