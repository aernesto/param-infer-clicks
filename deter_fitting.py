"""
The aim of this script is to infer the discounting & hazard rate used to produce simulated data.
"""
import pickle
# import copy
# import time
# import sys
from sympy import *
from official_fcns import *


def deter_fit(p):
    """
    computes the MSE and the average total admissible width for a given block size and block number.
    If the "trial_number" parameter (i.e. the block size) is 50 and the block number is 100,
    then the function will loop over a 100 blocks of size 50 and compute, across blocks, the MSE and
    the average total admissible width.
    :param p: a dict of parameters with the following keys:
        ['S','low_rate','high_rate','hazard_rate','T','best_gamma',
        'filename','samples_params','tot_trials_db','block_number',
        'trial_number','model_to_fit','reference_model','group_name']
    :return: (average scalar error across blocks, list of failure dicts)
    """
    with h5py.File(p['filename'], 'r') as db_file:
        group_name = p['group_name']
        # dset_trials = db_file[group_name + '/trials']
        # dset_info = db_file[group_name + '/trial_info']
        dset_dec_to_fit = db_file[group_name + '/decision_' + p['model_to_fit']]
        dset_dec_ref = db_file[group_name + '/decision_' + p['reference_model']]
        # print('\n-------------------\n{}{}'.format(p['model_to_fit'], p['reference_model']))
        T = p['T']
        if p['model_to_fit'] == 'lin':
            true_param = p['best_gamma']
        elif p['model_to_fit'] == 'nonlin':
            true_param = p['hazard_rate']
        else:
            raise ValueError("'model_to_fit' key in p argument has wrong value")

        running_scalar_error = 0
        bn = p['block_number']
        ntrials = p['trial_number']

        list_of_failure_dicts = []

        # inline function definition
        def get_block():
            indcs = np.random.choice(p['tot_trials_db'], size=ntrials, replace=False)  # sampling from bank
            doubled_list = list(zip(indcs, np.arange(ntrials), np.zeros(ntrials)))
            return np.array(doubled_list, dtype=[('index', 'i4'), ('order', 'i4'), ('new_order', 'i4')])
            # return slice(block_idx * ntrials, (block_idx + 1) * ntrials)

        for i in range(bn):
            block_slice = np.sort(get_block(), order='index')  # sort is because h5py requires increasing order
            reference_dec = dset_dec_ref[block_slice['index'], 0]
            decision_data = dset_dec_to_fit[block_slice['index'], 1:]

            # put trials back in order from sampling
            block_slice['new_order'] = np.arange(ntrials)
            retrieved_indices = np.sort(block_slice, order='order')
            reference_dec = reference_dec[retrieved_indices['new_order']]
            decision_data = decision_data[retrieved_indices['new_order']]

            all_sample_values, sp_tol = build_sample_vec(p['samples_params'])

            curr_err, failure_dict = get_block_err(reference_dec, decision_data, all_sample_values, sp_tol,
                                                   (p['samples_params']['start'], p['samples_params']['end']),
                                                   true_param)
            if curr_err is not None:
                running_scalar_error += curr_err / bn
            list_of_failure_dicts.append(failure_dict)
    return running_scalar_error, list_of_failure_dicts


def get_block_err(ref_dec, synthetic_dec, sample_array, sample_tolerance, sample_edges, true_parameter):
    """
    for a given block, computes the end admissible intervals, and from this list of intervals, computes
    the single scalar error defined as follows:
    1/ let f be a prob. density function that is uniform on the union intervals
    2/ let t be the true parameter and x live in the union of intervals
    3/ scalar_error = \int_{intervals} f(x)*(x-t)^2 dx
    :param ref_dec: value of reference decision for each trial
    :param synthetic_dec: matrix of decision per sample value
    :param sample_array: array of samples for parameter to fit
    :param sample_tolerance: width to add around valid samples
    :param sample_edges: tuple with lowest and highest sample values
    :param true_parameter: true parameter to recover
    :return: 2-tuple with scalar error and failure dict
    """
    fail_dict = {}

    # construct boolean matrix of compatible decisions
    col_vec_ref_dec = np.reshape(ref_dec, (-1, 1))
    bool_vec = col_vec_ref_dec == synthetic_dec

    # set places where synthetic decision was 0 to 1 in boolean matrix
    bool_vec[synthetic_dec == 0] = True

    # perform 'and' operation column-wise to get end-block valid samples
    # a valid sample is one that is compatible with all the trials in the block
    valid_samples = np.prod(bool_vec, axis=0)
    nvalid = np.sum(valid_samples)  # number of valid samples
    fail_dict['num_valid_samples'] = nvalid

    if nvalid:  # at least one sample value is compatible with every trial
        fail_dict['max_compatible_trials'] = bool_vec.shape[0]
        # print('{} trials case'.format(bool_vec.shape[0]))
        # print('{} valid samples found'.format(nvalid))
        # print('{} shape of sample array'.format(sample_array.shape))
        # print('{}, {} first and last original samples'.format(sample_array[0], sample_array[-1]))
        # print('{} sample edges'.format(sample_edges))
        sample_array[np.logical_not(valid_samples)] = np.nan
        # print('{} nb of nan samples'.format(np.sum(np.isnan(sample_array))))
    else:  # no sample is compatible with all trials
        interm = np.sum(bool_vec, axis=0)  # this is to avoid depletion
        nrows = np.max(interm)
        fail_dict['max_compatible_trials'] = nrows
        if nrows == 0:  # no sample is compatible with any trial
            sample_array[:] = np.nan
        else:  # some samples are valid for some, but not all, trials
            semi_valid_samples = interm == nrows
            fail_dict['num_semi_valid_samples'] = np.sum(semi_valid_samples)
            sample_array[np.logical_not(semi_valid_samples)] = np.nan

    interval_dict = {'interval': sample_edges, 'samples': sample_array}
    interval_list = reconstruct_interval(interval_dict, sample_tolerance)

    # fail_dict['final_samples'] = sample_array[np.logical_not(np.isnan(sample_array))]  # heavy, only for debug
    fail_dict['final_intervals'] = interval_list

    def lost_true_param():
        lost_bool = True
        for ivvl in interval_list:
            if ivvl[0] <= true_parameter <= ivvl[1]:
                lost_bool = False
                break
        return lost_bool

    fail_dict['lost_true_param'] = lost_true_param()
    return get_scalar_error_from_intervs(interval_list, true_parameter), fail_dict


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
    tot_trials = 10000
    block_number = 500
    file_list = [{'fname': 'data/S3lr5h1T2tr10000sp1000.h5', 'gname': 'lr5hr20h1T2', 'S': 3, 'lr': 5}]
    sple_dict = {'nonlin': {'start': 0, 'end': 10, 'number': 1000}, 'lin': {'start': 0, 'end': 10, 'number': 1000}}
    trial_report_list = [50, 100, 150, 200, 250, 300]
    params = {'hazard_rate': 1,
              'T': 2,
              'tot_trials_db': tot_trials}  # 100000}  # todo: read this off the db
    results = []
    for file in file_list:
        params['filename'] = file['fname']
        params['group_name'] = file['gname']
        params['S'] = file['S']
        params['low_rate'] = file['lr']
        params['high_rate'] = get_lambda_high(params['low_rate'], params['S'])
        if params['S'] in np.arange(0.5, 10.1, 0.5) and params['hazard_rate'] == 1:
            pol = False
        else:
            pol = True
        params['best_gamma'] = get_best_gamma(params['S'], params['hazard_rate'], polyfit=pol)

        report_values = {'linlin': [], 'nonlinnonlin': [], 'linnonlin': [], 'nonlinlin': []}
        for trial_report in trial_report_list:
            # print('trial {}'.format(trial_report))
            params['block_number'] = block_number# tot_trials // trial_report
            params['trial_number'] = trial_report
            for model_pair in [('lin', 'lin'), ('nonlin', 'nonlin'), ('lin', 'nonlin'), ('nonlin', 'lin')]:
                # print(''.join(model_pair))
                params['model_to_fit'], params['reference_model'] = model_pair
                params['samples_params'] = sple_dict[model_pair[0]]
                scalar_error, failures = deter_fit(params)
                # print(mse, avgwidth)
                report_values[''.join(model_pair)].append((scalar_error, failures))
                # print(report_values[''.join(model_pair)])
        results.append({'file': (file, trial_report_list), 'stats': report_values})
    # pickle.dump(results, open('/home/adrian/tosubmit_home/mse.pkl', 'wb'))
    pickle.dump(results, open('data/mse_local_S3lr5.pkl', 'wb'))
