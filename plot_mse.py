import pickle
import matplotlib
import matplotlib.pyplot as plt
font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 12}
lw = 4  # linewidth
matplotlib.rc('font', **font)

results = pickle.load(open("data/mse_local_new_error.pkl", "rb"))

f, ax_array = plt.subplots(4, 1, sharex='all',  # sharey='col',
                           squeeze=False)

# first file
fnum = 0
file = results[fnum]['file']
trial_numbers = file[1]

row2modl_map = ['linlin', 'nonlinnonlin', 'linnonlin', 'nonlinlin']
# loop over rows of suplot
for r in range(4):
    # lin 2 lin fit
    data = results[fnum]['stats'][row2modl_map[r]]

    # MSE plot
    ax = ax_array[r, 0]
    ax.plot(trial_numbers, data, linewidth=lw)
    ax.set_xlabel('trial nb')
    ax.set_ylabel('scalar error')
    ax.set_title(row2modl_map[r])

plt.tight_layout()
plt.show()
