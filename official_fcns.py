import numpy as np
import h5py


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


def get_best_gamma(s, h, polyfit=True):
    if polyfit:
        # coefficients of polynomial in decreasing power order (coefs[0]*x**n+...)
        coefs = [0.308812235133288,
                 0.781956794258376,
                 1.545682894736835]
        polyfcn = np.poly1d(coefs)
        return h*polyfcn(s/np.sqrt(h))
    else:
        corr = {'gamma': np.array([2.0848, 2.5828, 3.3143, 4.2789, 5.4162, 6.7457, 8.1371, 9.7199,
                                   11.3937, 13.2381, 15.1327, 17.2771, 19.5909, 22.0435, 24.6947,
                                   27.7241, 30.5711, 33.5354, 36.7966, 40.3143]),
                'S': np.arange(0.5, 10.1, 0.5)}
        iddx = np.where(corr['S'] == s)[0][0]
        return corr['gamma'][iddx]


def build_group_name(four_p):
    """
    :param four_p: (low click rate, high click rate, hazard rate, interrogation time)
    :return: string
    """
    #
    lowrate = ('{:f}'.format(round(four_p[0], 2))).rstrip('0').rstrip('.')
    highrate = ('{:f}'.format(round(four_p[1], 2))).rstrip('0').rstrip('.')
    hazard_rate = ('{:f}'.format(four_p[2])).rstrip('0').rstrip('.')
    interr_time = ('{:f}'.format(four_p[3])).rstrip('0').rstrip('.')
    return 'lr{}hr{}h{}T{}'.format(lowrate, highrate, hazard_rate, interr_time)


def get_lambda_high(lamb_low, s):
    """
    returns the high click rate, given the low click rate and S as input.
    :param lamb_low: low click rate
    :param s: S=(lambda_high-lambda_low)/sqrt(lambda_low+lambda_high)
    :return: value of lambda_high that fits with S and lambda_low
    """
    return (2 * lamb_low + s ** 2 + s * np.sqrt(s ** 2 + 8 * lamb_low)) / 2


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


def update_linear_decision_data(file_name, group_name, num_samples, sample_range, create_nonlin_db=False):
    """
    :param file_name: file name (string)
    :param group_name: group object from h5py module
    :param num_samples:
    :param sample_range: (starting value, ending value)
    :param create_nonlin_db:
    :return:
    """
    f = h5py.File(file_name, 'r+')
    group = f[group_name]
    info_dset = group['trial_info']
    trials_dset = group['trials']
    num_trials = trials_dset.shape[0]
    row_indices = range(num_trials)
    dset_name = 'decision_lin'
    if create_nonlin_db:
        # create dataset for nonlinear decisions
        group.create_dataset('decision_nonlin', (num_trials, num_samples),
                             dtype='i', maxshape=(100000, 10001))
    dset = group[dset_name]

    # store best gamma as attribute for future reference if doesn't exist
    skellam = info_dset.attrs['S']
    h = info_dset.attrs['h']
    if 'best_gamma' in dset.attrs.keys():
        best_gamma = dset.attrs['best_gamma']
    else:
        if skellam in np.arange(0.5, 10.1, 0.5) and h == 1:
            best_gamma = get_best_gamma(skellam, h, polyfit=False)
        else:
            best_gamma = get_best_gamma(skellam, h)
        dset.attrs['best_gamma'] = best_gamma
    gamma_samples, gamma_step = np.linspace(sample_range[0], sample_range[1], num_samples, retstep=True)
    attrslist = ['init_sample', 'end_sample', 'sample_step']
    values_dict = {'init_sample': sample_range[0],
                   'end_sample': sample_range[1],
                   'sample_step': gamma_step}
    for attrname in attrslist:
        if attrname not in dset.attrs.keys():
            dset.attrs[attrname] = values_dict[attrname]

    # populate dataset
    for row_idx in row_indices:
        stim = tuple(trials_dset[row_idx, :2])
        gamma_array = np.reshape(np.r_[best_gamma, gamma_samples], (-1, 1))
        dset[row_idx, :] = decide_linear(gamma_array, stim)
    f.flush()
    f.close()


def get_accuracy(filename, groupname, model, sample_col=0):
    """
    computes the accuracy of the given model in the given database, using the given sample value
    :param filename:
    :param groupname:
    :param model: either 'lin' or 'nonlin'
    :param sample_col:
    :return: accuracy between 0 and 1
    """
    if model not in ['lin', 'nonlin']:
        raise ValueError("model param should be either 'lin' or 'nonlin'")
    with h5py.File(filename, 'r') as f:
        grp = f[groupname]
        decision_dset = grp['decision_{}'.format(model)]
        info_dset = grp['trial_info']
        decisions = decision_dset[:, sample_col]
        correct_choice = info_dset[:, 1]  # end state of environment stored in 2nd column
        correct_decision = decisions == correct_choice
        return np.mean(correct_decision)


def reconstruct_interval(interval, tolerance):
    """
    splits an interval into sub-intervals
    :param interval: dict with: dict['interval'] = (lower_bound, upper_bound) and dict['samples'] = numpy.array
    :param tolerance: space between consecutive samples
    :return: list of sub-intervals, possibly empty
    """
    interval_list = []

    def appropriate_append(ivl):
        if not interval_list:
            interval_list.append(ivl)
        else:
            last_ivl = interval_list[-1]
            if last_ivl[1] == ivl[0]:
                interval_list[-1] = (last_ivl[0], ivl[1])
            else:
                interval_list.append(ivl)

    last_lower_bound, last_upper_bound = interval['interval']
    old_samples = interval['samples']
    lower_bound = None  # new lower bound
    last_up = 0 #last_lower_bound + tolerance  # new upper bound

    for indx, g in enumerate(old_samples):
        nan = np.isnan(g)
        if (lower_bound is None) and nan:
            continue
        elif g <= last_up:  # <= ?
            if g == last_lower_bound:
                lower_bound = last_lower_bound
            continue
        elif lower_bound is None:
            nlb = g - tolerance
            if nlb >= max(last_lower_bound, last_up):
                lower_bound = nlb
            else:
                lower_bound = max(last_lower_bound, last_up)
            if g == old_samples[-1]:  # if last sample reached
                appropriate_append((lower_bound, last_upper_bound))
            # for debug purposes
            if lower_bound < last_up:
                print('WARNING: overlapping consecutive intervals')
        elif nan:
            gm1 = old_samples[max(indx - 1, 0)]
            last_up = gm1 + tolerance
            appropriate_append((lower_bound, last_up))
            lower_bound = None
        # list of samples ended without nan values, so inherit upper bound from prev. setting
        elif g == old_samples[-1]:
            last_up = last_upper_bound
            appropriate_append((lower_bound, last_up))
    return interval_list


# todo: Fix function below
def get_scalar_error_from_intervs(list_of_intervals, true_param):
    if not list_of_intervals:
        # list is empty
        return None
    else:
        # tot_width = 0
        tot_intgl = 0
        for a, b in list_of_intervals:
            curr_width = b - a
            curr_intgl = (b**3 - a**3) / 3 + curr_width * true_param**2 + true_param * (a**2 - b**2)
            tot_intgl += curr_intgl
            # tot_width += curr_width
        try:
            scalar_err = tot_intgl  # / tot_width
        except ZeroDivisionError:
            print('Warning: samples depleted')
            return None
        return scalar_err


if __name__ == '__main__':
    # test for reconstruct interval
    # nan_cases = [[0,1,2,3,4],
    #              [0,1,2,3],
    #              [1,2,3,4],
    #              [0,1,3,4],
    #              [0,1,2,4],
    #              [0],
    #              [1],
    #              [3],
    #              [4],
    #              [0, 1],
    #              [0, 2],
    #              [0, 3],
    #              [0, 4],
    #              [1, 2],
    #              [1, 3],
    #              [1, 4],
    #              [2, 3],
    #              [2, 4],
    #              [3, 4],
    #              [0,1,2],
    #              [0,2,3],
    #              [0,2,4],
    #              [1,2,4],
    #              [1,2,2]]

    nan_cases = [range(4), range(8, 12), range(18, 22),range(28, 30)]
    indices = []
    for rr in nan_cases:
        for idx in rr:
            indices.append(idx)

    # print(indices)
    samples, tol = np.linspace(0, 10, 30, retstep=True)
    print('tolerance {:.2f}'.format(tol))
    samples[indices] = np.nan
    d={'interval': (0, 10), 'samples': samples}
    toprint = ['{:.2f}'.format(s) for s in samples]
    print('{}; {}'.format(toprint, reconstruct_interval(d, tol)))
