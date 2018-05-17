import pickle
import matplotlib
import matplotlib.pyplot as plt
font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 15}
lw = 4  # linewidth
matplotlib.rc('font', **font)

results = pickle.load(open("data/mse.pkl", "rb"))

f, ax_array = plt.subplots(2, 2, sharex='all', sharey='col', squeeze=False)

# first file
fnum = 0
file = results[fnum]['file']
trial_numbers = file[1]

row2modl_map = ['linnonlin', 'nonlinlin']
# loop over rows of suplot
for r in range(2):
    # lin 2 lin fit
    data = results[fnum]['stats'][row2modl_map[r]]

    # MSE plot
    ax = ax_array[r, 0]
    ax.plot(trial_numbers, [x[0] for x in data], linewidth=lw)
    ax.set_xlabel('trial nb')
    ax.set_ylabel('MSE')
    ax.set_title(row2modl_map[r])

    # Width plot
    ax = ax_array[r, 1]
    ax.plot(trial_numbers, [x[1] for x in data], linewidth=lw)
    ax.set_xlabel('trial nb')
    ax.set_ylabel('width')
    ax.set_title(row2modl_map[r])

plt.tight_layout()
plt.show()
