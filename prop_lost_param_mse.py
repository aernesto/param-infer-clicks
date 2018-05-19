import pickle
import matplotlib
import matplotlib.pyplot as plt
font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 12}
lw = 4  # linewidth
matplotlib.rc('font', **font)


def get_proportion_lost(list_of_dicts):
    plost = 0
    for d in list_of_dicts:
        if d['lost_true_param']:
            plost += 1
    return plost / len(list_of_dicts)


files = [{'filename': 'data/MSE_data_srvrs.pkl', 'num_files': 2, 'low_rate': [15, 1]},
         {'filename': 'data/mse_local_h20.pkl', 'num_files': 1, 'low_rate': [30]},
         {'filename': 'data/mse_new_local.pkl', 'num_files': 1, 'low_rate': [5]}]
row2modl_map = ['linlin', 'nonlinnonlin', 'linnonlin', 'nonlinlin']
proportions = {z: {} for z in row2modl_map}
for file_info in files:
    results = pickle.load(open(file_info['filename'], "rb"))
    for fnum in range(file_info['num_files']):
        lr = str(file_info['low_rate'][fnum])
        file = results[fnum]['file']
        trial_numbers = file[1]
        for r in row2modl_map:
            data = results[fnum]['stats'][r]
            proportions[r][lr] = [get_proportion_lost(d) for _, d in data]

# print proportions
for t in row2modl_map:
    for k in ['1', '5', '15', '30']:
        print('{} {} {}'.format(t, k, proportions[t][k]))
