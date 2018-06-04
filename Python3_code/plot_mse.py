import pickle
import matplotlib
import matplotlib.pyplot as plt
font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 12}
lw = 4  # linewidth
matplotlib.rc('font', **font)

results = pickle.load(open("data/mse_local_S3lr5.pkl", "rb"))

row2modl_map = ['linlin', 'nonlinnonlin', 'linnonlin', 'nonlinlin']


file = results[0]['file']
trial_numbers = file[1]
f, ax_array = plt.subplots(4, 1, sharex='all', sharey='col', squeeze=False)
# loop over rows of suplot
for r in range(4):
    # lin 2 lin fit
    data = results[0]['stats'][row2modl_map[r]]
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
