import sys
from official_fcns import *

import matplotlib.pyplot as plt
import pickle

if __name__ == '__main__':
    # If called from the command line with arguments
    if len(sys.argv) > 3:
        filename = sys.argv[1]
        grpname = sys.argv[2]
        model = sys.argv[3]
        try:
            splecol = int(sys.argv[4])
        except IndexError:
            splecol = 0
        print(get_accuracy(filename, grpname, model, splecol))
    else:
        fname = ["../data/S3lr5h1T2tr10000sp1000.h5",
                 "../data/S3lr2h1T2tr5Ksp10K.h5",
                 "/run/media/radillo/FreeAgent GoFlex Drive/data_clicks/srvr_data_1.h5"]
        gname = ["lr5hr20h1T2",
                 "lr2hr14h1T2",
                 "lr15hr36.5367250374h1T2"]
        ncols = [1000, 10000, 10000]
        models = ('lin', 'nonlin')

        accuracies = [np.zeros((2, ncols[i])) for i in range(3)]  # row1=nonlin; row0=lin
        for db in range(3):
            for col in np.arange(ncols[db]):
                accuracies[db][0, col] = get_accuracy(fname[db], gname[db], models[0], sample_col=col+1)
                accuracies[db][1, col] = get_accuracy(fname[db], gname[db], models[1], sample_col=col+1)

        # save data to file
        pickle.dump(accuracies, open('../data/accuracies_S3.pkl', 'wb'))

        # acc curves for linear model
        plt.figure()
        bg = get_best_gamma(3, 1, polyfit=False)
        plt.plot(np.linspace(0, 10, ncols[0]) / bg, accuracies[0][0, :],
                 np.linspace(0, 40, ncols[1]) / bg, accuracies[1][0, :],
                 np.linspace(0, 40, ncols[2]) / bg, accuracies[2][0, :],
                 linewidth=3)
        plt.plot([1, 1], [.6, .9])
        plt.title('lin accuracies')
        plt.ylabel('% correct')
        plt.xlabel('relative param')
        plt.legend(['S3lr5', 'S3lr2', 'S3lr15'])

        # acc curves for nonlinear model
        plt.figure()
        plt.plot(np.linspace(0, 10, ncols[0]), accuracies[0][1, :],
                 np.linspace(0, 40, ncols[1]), accuracies[1][1, :],
                 np.linspace(0, 40, ncols[2]), accuracies[2][1, :],
                 linewidth=3)
        plt.plot([1, 1], [.6, .9])
        plt.title('nonlin accuracies')
        plt.ylabel('% correct')
        plt.xlabel('relative param')
        plt.legend(['S3lr5', 'S3lr2', 'S3lr15'])

        plt.show()
