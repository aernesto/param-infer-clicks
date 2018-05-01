import pickle
import numpy as np
import matplotlib.pyplot as plt
f05 = pickle.load(open('data/fourBins_100000_lr15hr17.8664640249h1T2.pkl','rb'))
f8 = pickle.load(open('data/fourBins_100000_lr15hr101.258639865h1T2.pkl','rb'))
f3=pickle.load(open('data/fourBins_100000_lr15hr36.5367250374h1T2.pkl','rb'))
filelist = [f05, f3, f8]
for i in range(3):
    ff = filelist[i]
    means=np.array([np.mean(z) for z in ff])
    stdevs = np.array([np.std(z) for z in ff])
    plt.subplot(1, 3, i+1)
    plt.plot([.5,1,1.5,2], means, 'b')
    plt.plot([.5, 1, 1.5, 2], means-stdevs, 'r')
    plt.plot([.5, 1, 1.5, 2], np.minimum(means+stdevs, 40), 'r')
    plt.xlabel('last epoch duration')
    plt.ylabel('total width')
    plt.ylim(0, 40)
    plt.tight_layout()
plt.savefig('/home/radillo/Pictures/simulations/histograms/meanstdevs_2.svg', bbox_inches='tight')
