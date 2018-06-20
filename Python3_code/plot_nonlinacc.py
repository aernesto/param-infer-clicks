import pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 12}
# lw = 3  # linewidth
matplotlib.rc('font', **font)
ncols = [1000, 10000, 10000]
accuracies = pickle.load(open("../data/accuracies_S3.pkl", "rb"))

'''
the list accuracies contains accuracy data for the following files, in the following order
["../data/S3lr5h1T2tr10000sp1000.h5",
"../data/S3lr2h1T2tr5Ksp10K.h5",
"/run/media/radillo/FreeAgent GoFlex Drive/data_clicks/srvr_data_1.h5"]
'''

# for c in range(ncols[1]):
#     print('{:.6f}'.format(accuracies[1][0, c]))

'''
Alan's values
'''

assumed_h=np.arange(0.1,3.05,0.1)
max_acc_nonlin = [
   0.82889,
   0.84348,
   0.85142,
   0.85656,
   0.86000,
   0.86240,
   0.86406,
   0.86516,
   0.86591,
   0.86634,
   0.86651,
   0.86640,
   0.86613,
   0.86576,
   0.86525,
   0.86463,
   0.86395,
   0.86321,
   0.86238,
   0.86151,
   0.86060,
   0.85969,
   0.85871,
   0.85771,
   0.85665,
   0.85560,
   0.85453,
   0.85351,
   0.85243,
   0.85133]

# gamma_best=6.6449;
#
# assumed_gamma =[
#     0.38872
#     0.77745
#     1.16617
#     1.55490
#     1.94362
#     2.33235
#     2.72107
#     3.10979
#     3.49852
#     3.88724
#     4.27597
#     4.66469
#     5.05342
#     5.44214
#     5.83086
#     6.21959
#     6.60831
#     6.99704
#     7.38576
#     7.77449
#     8.16321
#     8.55193
#     8.94066
#     9.32938
#     9.71811
#    10.10683
#    10.49555
#    10.88428
#    11.27300
#    11.66173
#    12.05045
#    12.43918
#    12.82790
#    13.21662
#    13.60535
#    13.99407
#    14.38280
#    14.77152
#    15.16025
#    15.54897
#    15.93769
#    16.32642
#    16.71514
#    17.10387
#    17.49259
#    17.88132
#    18.27004
#    18.65876
#    19.04749
#    19.43621];
#
# max_acc_lin =[
#    0.74927
#    0.76570
#    0.78111
#    0.79536
#    0.80812
#    0.81925
#    0.82888
#    0.83693
#    0.84361
#    0.84903
#    0.85341
#    0.85673
#    0.85932
#    0.86111
#    0.86232
#    0.86306
#    0.86331
#    0.86322
#    0.86274
#    0.86205
#    0.86116
#    0.86007
#    0.85885
#    0.85749
#    0.85601
#    0.85445
#    0.85281
#    0.85105
#    0.84930
#    0.84752
#    0.84568
#    0.84384
#    0.84198
#    0.84013
#    0.83825
#    0.83632
#    0.83445
#    0.83260
#    0.83072
#    0.82886
#    0.82701
#    0.82517
#    0.82336
#    0.82155
#    0.81977
#    0.81795
#    0.81619
#    0.81445
#    0.81276
#    0.81107];

hvals=np.linspace(0, 40, ncols[2])
plt.plot(hvals[15:751], accuracies[2][1, 15:751],  # srvr_data_1
         assumed_h, max_acc_nonlin,
         linewidth=3)
plt.plot([1, 1], [.8, .9])
plt.title('nonlin accuracies; 100K trials')
plt.ylabel('% correct')
plt.xlabel('relative param')
plt.legend(['lr=15', 'lr=50'])

plt.show()
