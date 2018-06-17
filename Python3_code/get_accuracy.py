import sys
from official_fcns import *

import matplotlib.pyplot as plt

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
        fname = "../data/S3lr5h1T2tr10000sp1000.h5"
        gname = "lr5hr20h1T2"
        models = ('lin', 'nonlin')
        ncols = 1000
        accuracies = np.zeros((2, ncols))
        for col in np.arange(ncols):
            accuracies[0, col] = get_accuracy(fname, gname, models[0], sample_col=col+1)
            accuracies[1, col] = get_accuracy(fname, gname, models[1], sample_col=col+1)

        plt.plot(np.linspace(0, 10, ncols), accuracies[1, :],
                 np.linspace(0, 10, ncols)/get_best_gamma(3, 1, polyfit=False), accuracies[0, :],
                 linewidth=3)
        plt.plot([1, 1], [.65, .9])
        plt.ylabel('% correct')
        plt.xlabel('relative param')
        plt.legend(['nonlin', 'lin','true'])
        plt.show()
