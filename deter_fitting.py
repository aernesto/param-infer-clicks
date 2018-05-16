"""
The aim of this script is to infer the discounting & hazard rate used to produce simulated data.
"""

import matplotlib
# matplotlib.use('Agg')  # required on server to forbid X-windows usage
import pickle
import matplotlib.pyplot as plt
import copy
import time
import sys
import h5py
from sympy import *
from official_fcns import *

font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 15}

matplotlib.rc('font', **font)


# from fractions import Fraction
# from math import gcd

def deter_fit(p):
    """

    :param p: a dict of parameters with the following keys:
        ['S','low_rate','high_rate','hazard_rate','T','best_gamma',
        'filename','samples_params','tot_trials_db','block_number',
        'trial_number','model_to_fit','reference_model','group_name']
    :return: 2-tuple containing the average MSE across blocks and average total admissible width across blocks
    """
    with h5py.File(p['filename'], 'r') as db_file:
        group_name = p['group_name']
        # dset_trials = db_file[group_name + '/trials']
        # dset_info = db_file[group_name + '/trial_info']
        dset_dec_to_fit = db_file[group_name + '/decision_' + p['model_to_fit']]
        dset_dec_ref = db_file[group_name + '/decision_' + p['reference_model']]
        T = p['T']
        all_sample_values, sp_tol = build_sample_vec(p['samples_params'])
        if p['model_to_fit'] == 'lin':
            true_param = p['best_gamma']
        elif p['model_to_fit'] == 'nonlin':
            true_param = p['hazard_rate']
        else:
            raise ValueError("'model_to_fit' key in p argument has wrong value")

        running_mse = 0
        running_avg_width = 0
        N = p['block_number']
        ntrials = p['trial_number']

        # inline function definition
        def get_block(block_idx):
            return slice(block_idx * ntrials, (block_idx + 1) * ntrials)

        for i in range(N):
            block_slice = get_block(i)
            reference_dec = dset_dec_ref[block_slice, 0]
            decision_data = dset_dec_to_fit[block_slice, 1:]
            curr_sq_err, curr_width = get_block_width(reference_dec, decision_data, all_sample_values, sp_tol,
                                                      (p['samples_params']['start'], p['samples_params']['end']),
                                                      true_param)
            running_mse += curr_sq_err / N
            running_avg_width += curr_width / N

    return running_mse, running_avg_width


def get_block_width(ref_dec, synthetic_dec, sample_array, sample_tolerance, sample_edges, true_parameter):
    """
    :param ref_dec: value of reference decision for each trial
    :param synthetic_dec: matrix of decision per sample value
    :param sample_tolerance: width to add around valid samples
    :param sample_edges: tuple with lowest and highest sample values
    :param true_parameter: true parameter to recover
    :return: 2-tuple with squared error and total admissible width
    """
    # construct boolean matrix of compatible decisions
    col_vec_ref_dec = np.reshape(ref_dec, (-1, 1))
    bool_vec = col_vec_ref_dec == synthetic_dec

    # set places where synthetic decision was 0 to 1 in boolean matrix
    bool_vec[synthetic_dec == 0] = True

    # perform 'and' operation column-wise to get end-block valid samples
    valid_samples = np.prod(bool_vec, axis=0)

    # add tolerances around edges
    sample_array[np.logical_not(valid_samples)] = np.nan
    interval_dict = {'interval': sample_edges, 'samples': sample_array}
    interval_list = reconstruct_interval(interval_dict, sample_tolerance)

    # compute squared error and total width
    def get_sqerr_width_from_intervs(list_of_intervals):
        tot_width = 0
        tot_intgl = 0
        for ivl in list_of_intervals:
            curr_width = ivl[1] - ivl[0]
            curr_intgl = (ivl[1]**2-ivl[0]**2) / 2
            tot_intgl += curr_intgl
            tot_width += curr_width
        try:
            estimate = tot_intgl / tot_width
        except ZeroDivisionError:
            print('Warning: samples depleted')
            estimate = 0
        sqerr = (estimate - true_parameter)**2
        return sqerr, tot_width

    squared_error, total_width = get_sqerr_width_from_intervs(interval_list)

    return squared_error, total_width


def build_sample_vec(samples_params_dict):
    start = samples_params_dict['start']
    end = samples_params_dict['end']
    nb = samples_params_dict['number']
    return np.linspace(start, end, nb, retstep=True)


def get_best_gamma(skellam, h, polyfit=False):
    if polyfit:
        snr = skellam / np.sqrt(h)
        # todo: correct the flawed polynomial below
        return 1.45333 + 0.670241 * snr + 0.34324 * (snr ** 2) - 0.00275835 * (snr ** 3)
    else:
        corr = {'gamma': np.array([2.0848, 2.5828, 3.3143, 4.2789, 5.4162, 6.7457, 8.1371, 9.7199,
                                   11.3937, 13.2381, 15.1327, 17.2771, 19.5909, 22.0435, 24.6947,
                                   27.7241, 30.5711, 33.5354, 36.7966, 40.3143]),
                'S/sqrt(h)': np.arange(0.5, 10.1, 0.5)}
        iddx = np.where(corr['S/sqrt(h)'] == skellam / np.sqrt(h))[0][0]
        return corr['gamma'][iddx]


def dump_info(four_parameters, s, nt, nruns):
    print('S value: {}'.format(s))
    print('low click rate: {}'.format(four_parameters[0]))
    print('high click rate: {}'.format(four_parameters[1]))
    print('hazard rate: {}'.format(four_parameters[2]))
    print('interr. time: {}'.format(four_parameters[3]))
    print('nb of trials / hist: {}'.format(nruns))
    print('nb of trials in sequence: {}'.format(nt))


if __name__ == '__main__':
    """
    1/ Define parameters for synthetic data to read
    2/ Compute four sequences of admissible intervals for each trial (1 per model condition)
    3/ 
    """

    # 1/ Define parameters for synthetic data to read
    params = {'S': 3,
              'low_rate': 2,
              'hazard_rate': 1,
              'T': 2,
              'filename': 'data/S3lr2h1T2tr5Ksp10K.h5',
              'samples_params': {'start': 0, 'end': 40, 'number': 10000},
              'tot_trials_db': 5000,
              'block_number': 25,
              'trial_number': 200,
              'model_to_fit': 'nonlin',
              'reference_model': 'nonlin'}
    params['high_rate'] = get_lambda_high(params['low_rate'], params['S'])
    if params['S'] in np.arange(0.5, 10.1, 0.5) and params['hazard_rate'] == 1:
        pol = False
    else:
        pol = True
    params['best_gamma'] = get_best_gamma(params['S'], params['hazard_rate'], polyfit=pol)
    params['group_name'] = build_group_name((params['low_rate'],
                                             params['high_rate'],
                                             params['hazard_rate'],
                                             params['T']))
    start_time = time.time()
    MSE, AvgWidth = deter_fit(params)
    print('--{} seconds'.format(time.time() - start_time))
    print(MSE, AvgWidth)

