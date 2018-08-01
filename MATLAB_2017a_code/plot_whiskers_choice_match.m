clear
a=load('../data/choice_match_100_2.mat');
b=load('../data/choice_match_200_2.mat');
c=load('../data/choice_match_300_2.mat');
d=load('../data/choice_match_400_2.mat');
e=load('../data/choice_match_500_2.mat');
lw=4;
ms=8;
fs=20;
boxplot(e.match)
set(findobj(gca,'type','line'),'linew',lw)
set(gca,'linew',lw/2)
ylim([.79,.92])
ylabel('% match')
ax=gca;ax.FontSize=fs;
%saveas(gcf, ['whisker_',num2str(TN),'.png'])