clear
fname='data/S3lr2h1T2tr5Ksp10K.h5';
gname='/lr2hr14h1T2';
ntrials=5000;
nsamp=10000;
range_g=[0,40];
range_h=[0,40];
[l2l,l2n,n2l,n2n]=count_perfect_samples(fname,gname,ntrials,nsamp,...
    range_g(1),range_g(2),range_h(1), range_h(2));

ms=10;lw=1.5;fs=18;
xlim_lin=[min([min(l2l.values),min(l2n.values),l2l.true_param_fit]),...
    max([max(l2l.values),max(l2n.values),l2l.true_param_fit])];
xlim_nonlin=[min([min(n2l.values),min(n2n.values),n2l.true_param_fit]),...
    max([max(n2l.values),max(n2n.values),n2l.true_param_fit])];

subplot(2,2,1)
yval=l2l.maxprob*100;
plot(l2l.values,yval*ones(1,length(l2l.values)),'o','MarkerSize',ms,'LineWidth',lw)
hold on
plot(l2l.true_param_fit,yval,'*r','MarkerSize',ms,'LineWidth',lw)
hold off
ylabel('% trial')
ylim([yval-2,yval+2])
xlim(xlim_lin)
legend('sample','true','Location','southeast')
text(diff(xlim_lin)/2+xlim_lin(1),yval+1.5,[num2str(l2l.count),' values'],'FontSize',fs)
ax=gca;ax.FontSize=fs;
title('lin to lin')

subplot(2,2,2)
yval=l2n.maxprob*100;
plot(l2n.values,yval*ones(1,length(l2n.values)),'o','MarkerSize',ms,'LineWidth',lw)
hold on
plot(l2n.true_param_fit,yval,'*r','MarkerSize',ms,'LineWidth',lw)
hold off
ylim([yval-2,yval+2])
xlim(xlim_lin)
title('lin to nonlin')
legend('sample','true','Location','southeast')
text(diff(xlim_lin)/2+xlim_lin(1),yval+1.5,[num2str(l2n.count),' values'],'FontSize',fs)
ax=gca;ax.FontSize=fs;

subplot(2,2,3)
yval=n2l.maxprob*100;
plot(n2l.values,yval*ones(1,length(n2l.values)),'o','MarkerSize',ms,'LineWidth',lw)
hold on
plot(n2l.true_param_fit,yval,'*r','MarkerSize',ms,'LineWidth',lw)
hold off
ylabel('% trial')
ylim([yval-2,yval+2])
xlim(xlim_nonlin)
xlabel('param space')
title('nonlin to lin')
legend('sample','true','Location','southeast')
text(diff(xlim_nonlin)/2+xlim_nonlin(1),yval+1.5,[num2str(n2l.count),' values'],'FontSize',fs)
ax=gca;ax.FontSize=fs;

subplot(2,2,4)
yval=n2n.maxprob*100;
plot(n2n.values,yval*ones(1,length(n2n.values)),'o','MarkerSize',ms,'LineWidth',lw)
hold on
plot(n2n.true_param_fit,yval,'*r','MarkerSize',ms,'LineWidth',lw)
hold off
xlabel('param space')
ylim([yval-2,yval+2])
xlim(xlim_nonlin)
legend('sample','true','Location','southeast')
text(diff(xlim_nonlin)/2+xlim_nonlin(1),yval+1.5,[num2str(n2n.count),' values'],'FontSize',fs)
title('nonlin to nonlin')
ax=gca; ax.FontSize=fs;