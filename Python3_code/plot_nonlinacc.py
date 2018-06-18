import pickle
import numpy as np
import matplotlib.pyplot as plt
# font = {'family': 'DejaVu Sans',
#         'weight': 'bold',
#         'size': 12}
lw = 3  # linewidth
# matplotlib.rc('font', **font)
ncols = [1000, 10000, 10000]
accuracies = pickle.load(open("../data/accuracies_S3.pkl", "rb"))

for c in range(ncols[1]):
    print('{:.6f}'.format(accuracies[1][0, c]))


# plt.plot(np.linspace(0, 10, ncols[0]), accuracies[0][1, :],  # S3lr5
#          np.linspace(0, 40, ncols[1]), accuracies[1][1, :],  # S3lr2
#          np.linspace(0, 40, ncols[2]), accuracies[2][1, :],  # srvr_data_1
#          linewidth=3)
# plt.plot([1, 1], [.6, .9])
# plt.title('nonlin accuracies')
# plt.ylabel('% correct')
# plt.xlabel('relative param')
# plt.legend(['S3lr5', 'S3lr2', 'S3lr15'])
#
# plt.show()
