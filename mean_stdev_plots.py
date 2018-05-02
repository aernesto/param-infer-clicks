import pickle
import numpy as np
import matplotlib.pyplot as plt
files = [[] for _ in range(6)]
files[0] = (pickle.load(open('data/short_time_bins_lr15hr101.258639865h1T2.pkl', 'rb')), 8, 15)
files[1] = (pickle.load(open('data/short_time_bins_lr15hr17.8664640249h1T2.pkl', 'rb')), .5, 15)
files[2] = (pickle.load(open('data/short_time_bins_lr15hr36.5367250374h1T2.pkl', 'rb')), 3, 15)
files[3] = (pickle.load(open('data/short_time_bins_lr1hr11.6846584384h1T2.pkl', 'rb')), 3, 1)
files[4] = (pickle.load(open('data/short_time_bins_lr1hr1.84307033082h1T2.pkl', 'rb')), .5, 1)
files[5] = (pickle.load(open('data/short_time_bins_lr1hr66.941125497h1T2.pkl', 'rb')), 8, 1)

bin_edges = np.linspace(0, .8, 30+1)  # includes 0
bin_redges = bin_edges[1:]
for i in range(6):
    ff, S, lr = files[i]
    means = np.array([np.mean(z) for z in ff])
    stdevs = np.array([np.std(z) for z in ff])
    # if i == 0:
    #     plt.subplot(6, 1, 1)
    # elif i == 1:
    #     plt.subplot(6, 1, 3)
    # elif i == 2:
    #     plt.subplot(6, 1, 5)
    # elif i == 3:
    #     plt.subplot(6, 1, 2)
    # elif i == 4:
    #     plt.subplot(6, 1, 4)
    # elif i == 5:
    #     plt.subplot(6, 1, 6)
    plt.figure()
    plt.plot(bin_redges, means, 'b')
    plt.plot(bin_redges, means-stdevs, 'r')
    plt.plot(bin_redges, np.minimum(means+stdevs, 40), 'r')
    plt.title('S={}; low_rate={}'.format(S, lr))
    # if i in [2,5]:
    plt.xlabel('last epoch duration')
    # if i < 3:
    plt.ylabel('total width')
    plt.ylim(0, 40)
    plt.tight_layout()
    plt.savefig('/home/radillo/Pictures/simulations/whiskers/new_short_time_bins_S{}_lr{}.svg'.format(S, lr),
                bbox_inches='tight')
    plt.close()
