% plot and save whiskers of % match for each model pair

clear
a=load('../data/choice_match_100_4.mat');
b=load('../data/choice_match_200_4.mat');
c=load('../data/choice_match_300_4.mat');
d=load('../data/choice_match_400_4.mat');
e=load('../data/choice_match_500_4.mat');
files=[a,b,c,d,e];
lw=4;
ms=8;
fs=20;
nf=length(files);
TN=100:100:500;
for f=1:nf
    boxplot(files(f).match)
    set(findobj(gca,'type','line'),'linew',lw)
    set(gca,'linew',lw/2)
    hold on
    ax=gca;
    plot([ax.XLim(1), ax.XLim(2)],[.8787,.8787],'LineWidth',lw/2)
    plot([ax.XLim(1),ax.XLim(2)],[.86277248,0.86277248],'LineWidth',lw/2)
    hold off
    %ax.XTickLabels=['L-L','L-NL','NL-NL','NL-L'];
    ylim([.79,.92])
    ax.FontSize=fs;
    saveas(gcf, ['../figs/whisker_choice_',num2str(TN(f)),'.pdf'])
end