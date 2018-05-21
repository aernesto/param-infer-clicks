% testing time taken for pde solver on toy trial
% left_train = .1; % single click at .1 sec
% right_train = [];
% T = .5;  % trial lasts .5 sec
function pde_test
kappa = log(25/10); 
m = 0;
x = linspace(-2*kappa,2*kappa,40);
t = linspace(.1,.5,10);

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);

% A surface plot is often a good way to study a solution.
surf(x,t,u) 
title('Numerical solution computed with 20 mesh points.')
xlabel('Distance x')
ylabel('Time t')

% A solution profile can also be illuminating.
figure
plot(x,u(end,:))
title(['Solution at t = ',num2str(t(end))])
xlabel('Distance x')
ylabel(['u(x,',num2str(t(end)),')'])
% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
h=1;
c = 1;
f = 0;
s = 2*h*(cosh(x)*u+sinh(x)*DuDx);
% --------------------------------------------------------------
function u0 = pdex1ic(x)
kappa = log(25/10);
nsd = 1;
u0 = normpdf(x,kappa, nsd);
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = 0;
qr = 1;
