% plot accuracies
clear
noise01=load('../data/acc_noise_0.1.mat');
noise1=load('../data/acc_noise_1.mat');
noise2=load('../data/acc_noise_2.mat');
figure()
subplot(2,3,1)
plot(noise01.gammas,smooth(noise01.Acc))
subplot(2,3,4)
plot(noise01.h,noise01.Acc_h)
subplot(2,3,2)
plot(noise1.gammas,smooth(noise1.Acc))
subplot(2,3,5)
plot(noise1.h,noise1.Acc_h)
subplot(2,3,3)
plot(noise2.gammas,smooth(noise2.Acc))
subplot(2,3,6)
plot(noise2.h,noise2.Acc_h)