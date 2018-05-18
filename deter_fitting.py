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
        def get_block():
            indices = np.random.choice(p['tot_trials_db'], size=ntrials, replace=False)  # sampling from bank
            doubled_list = list(zip(indices, np.arange(ntrials), np.zeros(ntrials)))
            return np.array(doubled_list, dtype=[('index', 'i4'), ('order', 'i4'), ('new_order', 'i4')])
            # return slice(block_idx * ntrials, (block_idx + 1) * ntrials)

        for i in range(N):
            block_slice = np.sort(get_block(), order='index')  # sort is because h5py requires increasing order
            reference_dec = dset_dec_ref[block_slice['index'], 0]
            decision_data = dset_dec_to_fit[block_slice['index'], 1:]

            # put trials back in order from sampling
            block_slice['new_order'] = np.arange(ntrials)
            retrieved_indices = np.sort(block_slice, order='order')
            reference_dec = reference_dec[retrieved_indices['new_order']]
            decision_data = decision_data[retrieved_indices['new_order']]

            curr_sq_err, curr_width = get_block_width(reference_dec, decision_data, all_sample_values, sp_tol,
                                                      (p['samples_params']['start'], p['samples_params']['end']),
                                                      true_param)
            if (curr_sq_err is not None) and (curr_width is not None):
                running_mse += curr_sq_err / N
                running_avg_width += curr_width / N
        # if running_mse == 0 and running_avg_width == 0:
        #
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
    interm = np.cumprod(bool_vec, axis=0)  # this is to avoid depletion
    nonzero_row_indices = np.nonzero(interm)[0]
    if nonzero_row_indices.size > 0:
        valid_samples = interm[np.max(nonzero_row_indices), :]  # last nonzero row
        # add tolerances around edges
        sample_array[np.logical_not(valid_samples)] = np.nan
    else:
        valid_samples = np.array([])
        sample_array[:] = np.nan

    interval_dict = {'interval': sample_edges, 'samples': sample_array}
    interval_list = reconstruct_interval(interval_dict, sample_tolerance)

    # compute squared error and total width
    def get_sqerr_width_from_intervs(list_of_intervals):
        if not list_of_intervals:
            # list is empty
            return None, None
        else:
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

    return get_sqerr_width_from_intervs(interval_list)


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
    trial_report_list = [50, 100, 150, 200]
    params = {'hazard_rate': 1,
              'T': 2,
              'samples_params': {'start': 0, 'end': 10, 'number': 1000},
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

                mse, avgwidth = deter_fit(params)
                # print(mse, avgwidth)
                report_values[''.join(model_pair)].append((mse, avgwidth))
                # print(report_values[''.join(model_pair)])
        results.append({'file': (file, trial_report_list), 'stats': report_values})
    # pickle.dump(results, open('/home/adrian/tosubmit_home/mse.pkl', 'wb'))
    pickle.dump(results, open('data/mse_local_sampling.pkl', 'wb'))
