"""
The aim of this script is to infer the discounting & hazard rate used to produce simulated data.
"""
import pickle
# import copy
# import time
# import sys
from sympy import *
from official_fcns import *


def get_intervs(p):
    """
    computes the admissible intervals for a given block of trials.
    :param p: a dict of parameters with the following keys:
        ['S','low_rate','high_rate','hazard_rate','T','best_gamma',
        'filename','samples_params','tot_trials_db','block_number',
        'trial_number','model_to_fit','reference_model','group_name']
    :return: (lower bound, upper bound)
    """
    with h5py.File(p['filename'], 'r') as db_file:
        group_name = p['group_name']
        dset_dec_to_fit = db_file[group_name + '/decision_' + p['model_to_fit']]
        dset_dec_ref = db_file[group_name + '/decision_' + p['reference_model']]
        T = p['T']
        if p['model_to_fit'] == 'lin':
            true_param = p['best_gamma']
        elif p['model_to_fit'] == 'nonlin':
            true_param = p['hazard_rate']
        else:
            raise ValueError("'model_to_fit' key in p argument has wrong value")

        bn = p['block_number']
        ntrials = p['trial_number']

        for i in range(bn):
            reference_dec = dset_dec_ref[:ntrials, 0]
            decision_data = dset_dec_to_fit[:ntrials, 1:]

            all_sample_values, sp_tol = build_sample_vec(p['samples_params'])

            intervs = get_block_ivs(reference_dec, decision_data, all_sample_values, sp_tol,
                                    (p['samples_params']['start'], p['samples_params']['end']),
                                    true_param)
    return intervs


def get_block_ivs(ref_dec, synthetic_dec, sample_array, sample_tolerance, sample_edges, true_parameter):
    """
    for a given block, computes the end admissible intervals
    :param ref_dec: value of reference decision for each trial
    :param synthetic_dec: matrix of decision per sample value
    :param sample_array: array of samples for parameter to fit
    :param sample_tolerance: width to add around valid samples
    :param sample_edges: tuple with lowest and highest sample values
    :param true_parameter: true parameter to recover
    :return: 2-tuple with scalar error and failure dict
    """

    # construct boolean matrix of compatible decisions
    col_vec_ref_dec = np.reshape(ref_dec, (-1, 1))
    bool_vec = col_vec_ref_dec == synthetic_dec

    # set places where synthetic decision was 0 to 1 in boolean matrix
    bool_vec[synthetic_dec == 0] = True

    # perform 'and' operation column-wise to get end-block valid samples
    # a valid sample is one that is compatible with all the trials in the block
    valid_samples = np.prod(bool_vec, axis=0)
    nvalid = np.sum(valid_samples)  # number of valid samples

    if nvalid:  # at least one sample value is compatible with every trial
        sample_array[np.logical_not(valid_samples)] = np.nan
    else:  # no sample is compatible with all trials
        interm = np.sum(bool_vec, axis=0)  # this is to avoid depletion
        nrows = np.max(interm)
        if nrows == 0:  # no sample is compatible with any trial
            sample_array[:] = np.nan
        else:  # some samples are valid for some, but not all, trials
            semi_valid_samples = interm == nrows
            sample_array[np.logical_not(semi_valid_samples)] = np.nan

    interval_dict = {'interval': sample_edges, 'samples': sample_array}

    return reconstruct_interval(interval_dict, sample_tolerance)


def build_sample_vec(samples_params_dict):
    start = samples_params_dict['start']
    end = samples_params_dict['end']
    nb = samples_params_dict['number']
    return np.linspace(start, end, nb, retstep=True)


def dump_info(four_parameters, s, nt, nruns):
    print('S value: {}'.format(s))
    print('low click rate: {}'.format(four_parameters[0]))
    print('high click rate: {}'.format(four_parameters[1]))
    print('hazard rate: {}'.format(four_parameters[2]))
    print('interr. time: {}'.format(four_parameters[3]))
    print('nb of trials / hist: {}'.format(nruns))
    print('nb of trials in sequence: {}'.format(nt))


if __name__ == '__main__':
    tot_trials = 10000
    block_number = 1
    file_list = [{'fname': 'data/S3lr5h1T2tr10000sp1000.h5', 'gname': 'lr5hr20h1T2', 'S': 3, 'lr': 5}]
    sple_dict = {'nonlin': {'start': 0, 'end': 10, 'number': 1000}, 'lin': {'start': 0, 'end': 10, 'number': 1000}}
    trial_report_list = [50]
    params = {'filename': 'data/S3lr5h1T2tr10000sp1000.h5',
              'group_name': 'lr5hr20h1T2',
              'S': 3, 'low_rate': 5, 'high_rate': 20,
              'best_gamma': get_best_gamma(3, 1, polyfit=False),
              'hazard_rate': 1,
              'T': 2,
              'tot_trials_db': tot_trials,
              'block_number': block_number,
              'trial_number': 50}
    for model_pair in [('lin', 'lin'), ('nonlin', 'nonlin'), ('lin', 'nonlin'), ('nonlin', 'lin')]:
        # print(''.join(model_pair))
        params['model_to_fit'], params['reference_model'] = model_pair
        params['samples_params'] = sple_dict[model_pair[0]]
        print(model_pair, get_intervs(params))
