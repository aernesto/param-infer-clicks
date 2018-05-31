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


files = [{'filename': 'data/mse_local_S3lr5.pkl', 'num_files': 1, 'low_rate': [5]}]
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
    for k in ['5']:
        print('{} {} {}'.format(t, k, proportions[t][k]))
