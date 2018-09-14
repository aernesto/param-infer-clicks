clear
% best gamma data
bg=[2.0848, 2.5828, 3.3143, 4.2789, 5.4162, 6.7457, 8.1371, 9.7199,...
 11.3937, 13.2381, 15.1327, 17.2771, 19.5909, 22.0435, 24.6947,...
 27.7241, 30.5711, 33.5354, 36.7966, 40.3143];
% corresponding values of S/sqrt(h)
sh=0.5:0.5:10;
p2=polyfit(sh, bg, 2); y2=polyval(p2,sh);
p3=polyfit(sh, bg, 3); y3=polyval(p3,sh);
p4=polyfit(sh, bg, 4); y4=polyval(p4,sh);
plot(sh,bg,'o')
hold on
plot(sh,y2,sh,y3,sh,y4)
legend('data','2','3','4')
xlabel('S/sqrt(h)')
ylabel('best gamma')
ax = gca;
ax.FontSize=20;
hold off
