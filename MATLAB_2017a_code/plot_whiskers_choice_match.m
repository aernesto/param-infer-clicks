clear
a=load('../data/choice_match_100_1.mat');
b=load('../data/choice_match_200_1.mat');
c=load('../data/choice_match_300_1.mat');
d=load('../data/choice_match_400_1.mat');
e=load('../data/choice_match_500_1.mat');
lw=4;
ms=8;
fs=20;
boxplot(e.match)
set(findobj(gca,'type','line'),'linew',lw)
set(gca,'linew',lw/2)
%ylim([.85,.9])
ylabel('% match')
ax=gca;ax.FontSize=fs;
%saveas(gcf, ['whisker_',num2str(TN),'.png'])