"""
The aim of this script is to manage the databases in hdf5 format
"""
import numpy as np
# import copy
import time
import sys
import h5py

# from sympy import *

# from fractions import Fraction
# from math import gcd


"""
----------------------------GENERAL PURPOSE FUNCTIONS
"""


# def end_point_nonlin(init_cond, end_time, hr):
#     if init_cond == 0:
#         return 0
#     else:
#         return float(2*mpmath.acoth(np.exp(2*hr*end_time)*mpmath.coth(init_cond/2)))


def get_lambda_high(lamb_low, s):
    """
    returns the high click rate, given the low click rate and S as input.
    :param lamb_low: low click rate
    :param s: S=(lambda_high-lambda_low)/sqrt(lambda_low+lambda_high)
    :return: value of lambda_high that fits with S and lambda_low
    """
    return (2 * lamb_low + s ** 2 + s * np.sqrt(s ** 2 + 8 * lamb_low)) / 2


"""
----------------------------FUNCTIONS USED TO SIMULATE DATA
"""


def gen_cp(duration, rate):
    """
    generate the CP times from the trial by successively sampling
    from an exponential distribution
    :param duration: trial duration in sec
    :param rate: hazard rate in Hz
    :return: numpy array of change point times strictly positive and inferior to duration
    """
    cp_times = []
    cp_time = 0
    while cp_time < duration:
        dwell_time = np.random.exponential(1. / rate)
        cp_time += dwell_time
        if cp_time < duration:
            cp_times += [cp_time]

    return np.array(cp_times)


def gen_stim(ct, rate_l, rate_h, dur):
    """
    generates two simultaneous clicks streams (one per ear), selecting the
    initial environmental state at random (P_0(H) is uniform).
    Recall that environment == 1 means that high rate is presented to right ear
    -1 otherwise.
    :param ct: numpy array of change point times
    :param rate_l: low click rate in Hz
    :param rate_h: high click rate in Hz
    :param dur: trial duration in seconds
    :return: tuple with 2 elements,
     The first element is itself the tuple (left_stream, right_stream) (each stream is numpy array)
     The second element is the last environmental state
    """
    state = np.zeros(len(ct) + 1)
    num_trains = len(state)  # number of trains to stack, for each ear
    state[0] = np.random.choice([-1, 1])

    # flip state after each change point
    for counter in range(1, num_trains):
        state[counter] = -1 * state[counter - 1]

    left_stream = []  # storing click trains for each ear
    right_stream = []

    # construct trains between each change point
    for tt in range(num_trains):
        # extract time length of current train
        if tt == 0:
            if len(ct) > 0:
                time_length = ct[tt]
                offset = 0
            else:
                time_length = dur
                offset = 0
        elif tt == (num_trains - 1):
            offset = ct[-1]
            time_length = dur - offset
        else:
            offset = ct[tt - 1]
            time_length = ct[tt] - offset

        # construct trains for both ears, depending on envt state
        left_train_low = [XX + offset for XX in gen_cp(duration=time_length, rate=rate_l)]
        left_train_high = [XX + offset for XX in gen_cp(duration=time_length, rate=rate_h)]
        right_train_high = [XX + offset for XX in gen_cp(duration=time_length, rate=rate_h)]
        right_train_low = [XX + offset for XX in gen_cp(duration=time_length, rate=rate_l)]

        if state[tt] == 1:  # evaluates to true if envt is in state H+ ---> high rate to right ear
            left_stream += left_train_low
            right_stream += right_train_high
        else:  # envt in state H- ---> high rate to left ear
            left_stream += left_train_high
            right_stream += right_train_low

    return (np.array(left_stream), np.array(right_stream)), state[0], state[-1]


def decide_linear(gammas_array, stimulus, init_cond=0):
    """
    computes decision using the linear model
    :param gammas_array: numpy n-by-1 array of linear discounting terms
    :param stimulus: 2-tuple for left and right click trains
    :param init_cond: initial value of accumulation variable
    :return: -1, 0 or 1 for left, undecided and right
    """
    y = init_cond
    gammas = gammas_array.reshape((-1, 1))

    # right train
    right_train = stimulus[1].reshape((1, -1))
    y += np.exp(gammas @ right_train).sum(axis=1)

    # left train
    left_train = stimulus[0].reshape((1, -1))
    y -= np.exp(gammas @ left_train).sum(axis=1)

    return np.sign(y)


def create_hfd5_data_structure(file, groupname, num_trials):
    """
    As of now, decision datasets contain 10001 columns because first col contains
    decision for true parameter and the remaining 10,000 correspond to decision data
    parameter samples
    Current max nb of trials per dataset = 100,000
    :param file: h5py.File
    :param groupname:
    :param num_trials: nb of trials
    :return: created group
    """
    max_trials = 100000
    num_samples = 10000
    group = file.create_group(groupname)
    dt = h5py.special_dtype(vlen=np.dtype('f'))
    group.create_dataset('trials', (num_trials, 3), maxshape=(max_trials, 10), dtype=dt)
    group.create_dataset('trial_info', (num_trials, 3), maxshape=(max_trials, 10), dtype='f')
    group.create_dataset('decision_lin', (num_trials, num_samples+1), dtype='i', maxshape=(max_trials, num_samples+1))
    group.create_dataset('decision_nonlin', (num_trials, num_samples+1), dtype='i', maxshape=(max_trials, num_samples+1))
    return group


def build_group_name(four_p):
    """
    :param four_p: (low click rate, high click rate, hazard rate, interrogation time)
    :return: string
    """
    lowrate = str(round(four_p[0], 2))
    highrate = str(round(four_p[1], 2))
    hazard_rate = str(four_p[2])
    interr_time = str(four_p[3])
    return 'lr{}hr{}h{}T{}'.format(lowrate, highrate, hazard_rate, interr_time)


def populate_hfd5_db(fname, four_par, num_of_trials):
    """generate stimulus data and store as hdf5 file"""
    # open/create file
    f = h5py.File(fname, 'a')
    ll, lh, h, t = four_par
    # get/create group corresponding to parameters
    group_name = build_group_name((ll, lh, h, t))
    if group_name in f:  # if dataset already exists, only expand it with new data
        grp = f[group_name]

        # get datasets
        trials_data = grp['trials']
        info_data = grp['trial_info']

        # resizing operation before inserting new data:
        old_size = trials_data.len()
        new_size = old_size + num_of_trials
        trials_data.resize(new_size, axis=0)
        info_data.resize(new_size, axis=0)

        # get row indices of new data to insert
        row_indices = np.r_[old_size:new_size]

        # version number of new data to insert
        data_version = info_data.attrs['last_version'] + 1
    else:  # if dataset doesn't exist, create it
        grp = create_hfd5_data_structure(f, group_name, num_of_trials)

        # get trials dataset
        trials_data = grp['trials']
        # get row indices of new data to insert
        row_indices = np.r_[:num_of_trials]

        # create info on data
        info_data = grp['trial_info']  # info dataset
        info_data.attrs['h'] = h
        info_data.attrs['T'] = t
        info_data.attrs['low_click_rate'] = ll
        info_data.attrs['high_click_rate'] = lh
        info_data.attrs['S'] = (lh - ll) / np.sqrt(ll + lh)
        data_version = 1  # version number of new data to insert

    # populate database
    for row_idx in row_indices:
        # vector of CP times
        cptimes = gen_cp(t, h)
        trials_data[row_idx, 2] = cptimes

        # stimulus (left clicks, right clicks)
        (left_clicks, right_clicks), init_state, end_state = gen_stim(cptimes, ll, lh, t)
        trials_data[row_idx, :2] = left_clicks, right_clicks

        # populate info dataset
        info_data[row_idx, :] = init_state, end_state, data_version

    info_data.attrs['last_version'] = data_version
    f.flush()
    f.close()


def get_best_gamma(ratio_rates, h, polyfit=False):
    if polyfit:
        snr = ratio_rates / np.sqrt(h)
        # todo: correct the flawed polynomial below
        return 1.45333 + 0.670241 * snr + 0.34324 * (snr ** 2) - 0.00275835 * (snr ** 3)
    else:
        corr = {'gamma': np.array([2.0848, 2.5828, 3.3143, 4.2789, 5.4162, 6.7457, 8.1371, 9.7199,
                                   11.3937, 13.2381, 15.1327, 17.2771, 19.5909, 22.0435, 24.6947,
                                   27.7241, 30.5711, 33.5354, 36.7966, 40.3143]),
                'S/h': np.arange(0.5, 10.1, 0.5)}
        iddx = np.where(corr['S/h'] == ratio_rates)[0][0]
        return corr['gamma'][iddx]


def update_linear_decision_data(file_name, group_name, verbose=False):
    """
    :param file_name: file name (string)
    :param group_name: group object from h5py module
    :return:
    """
    f = h5py.File(file_name, 'a')
    group = f[group_name]
    info_dset = group['trial_info']
    trials_dset = group['trials']
    num_trials = trials_dset.shape[0]
    row_indices = range(num_trials)
    dset_name = 'decision_lin'
    num_samples = 10000
    # check whether decision dataset exists
    if dset_name in group:
        dset = group[dset_name]
        dset_length, dset_width = dset.shape
        # resize for new data
        if dset_length < num_trials:
            dset.resize(num_trials, axis=0)
            row_indices = np.arange(dset_length, num_trials)
        else:
            if verbose:
                print('possibly overwriting data')

        if dset_width < num_samples + 1:
            dset.resize(num_samples+1, axis=1)
        else:
            if verbose:
                print('possibly overwriting data')
    else:
        dset = group.create_dataset(dset_name, (num_trials, num_samples+1), dtype='i', maxshape=(100000, num_samples+1))
        row_indices = np.arange(num_trials)

    # store best gamma as attribute for future reference
    best_gamma = get_best_gamma(round(info_dset.attrs['S'], 4), 1)
    dset.attrs['best_gamma'] = best_gamma
    gamma_samples, gamma_step = np.linspace(0, 40, 10000, retstep=True)
    dset.attrs['init_sample'] = gamma_samples[0]
    dset.attrs['end_sample'] = gamma_samples[-1]
    dset.attrs['sample_step'] = gamma_step

    # populate dataset
    for row_idx in row_indices:
        # get trial object
        stim = tuple(trials_dset[row_idx, :2])
        gamma_array = np.reshape(np.r_[best_gamma,gamma_samples], (-1,1))
        dset[row_idx, :] = decide_linear(gamma_array, stim)
    f.flush()
    f.close()


def dump_info(four_parameters, s, nt, nruns):
    print('S value: {}'.format(s))
    print('low click rate: {}'.format(four_parameters[0]))
    print('high click rate: {}'.format(four_parameters[1]))
    print('hazard rate: {}'.format(four_parameters[2]))
    print('interr. time: {}'.format(four_parameters[3]))
    print('nb of trials / hist: {}'.format(nruns))
    print('nb of trials in sequence: {}'.format(nt))


if __name__ == '__main__':
    # parameter vectors
    a_S = np.arange(0.5, 6.1, 0.5)

    a_ll = (1, 15, 30)

    # scalar parameters
    int_time = 2
    hazard = 1
    filename = 'data/test_light.h5'
    start_time = time.time()
    number_of_trials = 100

    start_time = time.time()
    for S in a_S:
        for lambda_low in a_ll:  # low click rate
            lambda_high = get_lambda_high(lambda_low, S)
            four_params = (lambda_low, lambda_high, hazard, int_time)
            grp_name = build_group_name(four_params)
            true_g = get_best_gamma(S, hazard)

            populate_hfd5_db(filename, four_params, number_of_trials)
            update_linear_decision_data(filename, grp_name, verbose=True)

    print("--- {} seconds ---".format(time.time() - start_time))
