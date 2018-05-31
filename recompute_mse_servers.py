import pickle
from official_fcns import *
import matplotlib
import matplotlib.pyplot as plt
font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 12}
lw = 4  # linewidth
matplotlib.rc('font', **font)

results = pickle.load(open("data/MSE_data_srvrs.pkl", "rb"))

"""
results is a list of two dicts, each having same structure.
The fields of each dict are:
"file" and "stats"

The value corresponding to "file" is a 2-tuple.
the following is an example of such tuple:
({'fname': '/storage/adrian/srvr_data_1.h5', 'gname': 'lr15hr36.5367250374h1T2', 'S': 3, 'lr': 15}, 
[50, 100, 150, 200, 250, 300, 350, 400])

The value corresponding to "stats" is a dict with the following keys:
'linlin', 'nonlinnonlin', 'linnonlin', 'nonlinlin'
The value corresponding to each one of these keys is a list of 8 2-tuples, one per trial number.
Each 2-tuple looks like the following:
(wrong scalar error, list of N dicts where N is the number of trial blocks)

A single one of the dicts mentioned in the last line contains a key "final_intervals". 
The value for this key is a list of admissible intervals.


"""

row2modl_map = ['linlin', 'nonlinnonlin', 'linnonlin', 'nonlinlin']
fnum = 1  # 1 for data_S_2_5.h5 and 0 for srvr_data_1.h5
file = results[fnum]['file']
trial_numbers = file[1]
f, ax_array = plt.subplots(4, 1, sharex='all', sharey='col', squeeze=False)
# loop over rows of suplot
for r in range(4):
    # true parameter
    if row2modl_map[r] == 'linlin' or row2modl_map[r] == 'linnonlin':
        best_gamma = get_best_gamma(file[0]['S'], 1, polyfit=False)
        tp = best_gamma
    else:
        tp = 1
    # lin 2 lin fit
    data = results[fnum]['stats'][row2modl_map[r]]
    # explore failures
    # fail_dicts = [x[1] for x in data]
    # MSE plot
    ax = ax_array[r, 0]
    lldicts = [x[1] for x in data]  # each element is a list of 1000 dicts
    scalar_errors = [np.mean([get_scalar_error_from_intervs(x['final_intervals'], tp) for x in l]) for l in lldicts]
    ax.plot(trial_numbers, scalar_errors, 'o--', markersize=10, linewidth=lw)
    ax.set_xlabel('trial nb')
    ax.set_ylabel('scalar error')
    ax.set_title(row2modl_map[r])
plt.tight_layout()
plt.show()
