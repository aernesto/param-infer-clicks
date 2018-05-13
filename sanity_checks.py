from official_fcns import *
import matplotlib.pyplot as plt


if __name__ == "__main__":
    fname = []
    gname =[]
#    fname+=['data/S2lr3h1T2tr5000sp5.h5']; gname+=['lr3hr10.29h1T2'] 
#    fname+=['data/S2lr10h1T2tr5000sp5.h5']; gname+=['lr10hr21.17h1T2']
#    fname+=['data/S2lr25h1T2tr5000sp5.h5']; gname+=['lr25hr41.28h1T2']
    fname += ['data/S6lr25h1T2tr1000sp5.h5']; gname += ['lr25hr89.09h1T2']
    fname += ['data/S6lr10h1T2tr1000sp5.h5']; gname += ['lr10hr60.31h1T2']
    fname += ['data/S6lr3h1T2tr1000sp5.h5']; gname += ['lr3hr44.24h1T2']
    lin_acc=[]
    nonlin_acc=[]
    for i in range(len(fname)):
        lin_acc.append(float(get_accuracy(fname[i],gname[i],'lin')))
        nonlin_acc.append(float(get_accuracy(fname[i],gname[i],'nonlin')))
    print('lin accuracies')
    print(lin_acc)
    print('nonlin accuracies')
    print(nonlin_acc)
    plt.plot([3,10,25],lin_acc,'ob')
    plt.plot([3,10,25],nonlin_acc,'*r')
    plt.show()

