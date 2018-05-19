import pickle
import matplotlib
import matplotlib.pyplot as plt
font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 12}
lw = 4  # linewidth
matplotlib.rc('font', **font)

results = pickle.load(open("data/MSE_data_srvrs.pkl", "rb"))

row2modl_map = ['linlin', 'nonlinnonlin', 'linnonlin', 'nonlinlin']

for fnum in range(2):
    if fnum == 0:
        continue
    file = results[fnum]['file']
    trial_numbers = file[1]
    f, ax_array = plt.subplots(4, 1, sharex='all', sharey='col', squeeze=False)
    # loop over rows of suplot
    for r in range(4):
        # lin 2 lin fit
        data = results[fnum]['stats'][row2modl_map[r]]
        # explore failures
        # fail_dicts = [x[1] for x in data]

        # MSE plot
        ax = ax_array[r, 0]
        ax.plot(trial_numbers, [x[0] for x in data], 'o--', markersize=10, linewidth=lw)
        ax.set_xlabel('trial nb')
        ax.set_ylabel('scalar error')
        ax.set_title(row2modl_map[r])
    plt.tight_layout()
    plt.show()
