from official_fcns import *
#import matplotlib.pyplot as plt


if __name__ == "__main__":
    gname = ['lr15hr28.14h1T2','lr15hr57.6h1T2','lr1hr27.86h1T2','lr1hr6.46h1T2','lr30hr47.62h1T2','lr30hr83.2h1T2']
    fname = '/storage/adrian/data_S_2_5.h5' 
    lin_acc = []
    nonlin_acc = []
    for i in range(len(gname)):
        lin_acc.append(float(get_accuracy(fname,gname[i],'lin')))
        nonlin_acc.append(float(get_accuracy(fname,gname[i],'nonlin')))
    print('lin accuracies')
    print('S=2: {0[3]} {0[0]} {0[4]}; S=5: {0[2]} {0[1]} {0[5]}'.format(lin_acc))
    print('nonlin accuracies')
    print('S=2: {0[3]} {0[0]} {0[4]}; S=5: {0[2]} {0[1]} {0[5]}'.format(nonlin_acc))

