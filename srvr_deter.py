"""
The aim of this script is to infer the discounting rate used to produce simulated data.
The script contains two types of functions.
    The ones used to simulate the data.
    The ones used to infer the parameter.

"""
import numpy as np
import matplotlib
matplotlib.use('Agg')  # required on server to forbid X-windows usage
import matplotlib.pyplot as plt
import copy
import time
import sys
import h5py

font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 22}

# matplotlib.rc('font', **font)
# from fractions import Fraction
# from math import gcd

"""
----------------------------CLASS DEFINITIONS
"""


class Trial:
    def __init__(self, stimulus, gamma, trial_number, init_gammas=None, init_gamma_range=(0, 40),
                 tolerance=.05, verbose=False, hazard_rate=1, decision_datum=None):
        """
        :param stimulus: 2-tuple of click trains (left, right)
        :param gamma: true gamma with which decision data should be computed
        :param trial_number: number within external sequence
        :param init_gammas: initial intervals of admissible gammas (list of dicts)
        :param tolerance: max distance allowed between two consecutive samples
        """
        self.number = trial_number
        self.stimulus = stimulus
        self.true_gamma = gamma
        self.tolerance_gamma = tolerance
        if decision_datum is None:
            self.decision = self.decide(np.array([self.true_gamma]))  # -1 for left, 1 for right, 0 for undecided
        else:
            self.decision = decision_datum
        self.nonlin_dec = None  # todo: improve this
        if init_gammas is None:
            self.admissible_gammas = self.gen_init_gammas(init_gamma_range)
        else:
            self.admissible_gammas = init_gammas
        self.refined = False  # True if self.refine_admissible_gammas() has been called
        self.total_width = self.get_total_width()
        self.num_intervals = len(self.admissible_gammas)
        self.verbose = verbose
        if verbose:
            self.report('\n TRIAL {} INSTANTIATED \n'.format(self.number))

    def report(self, string='', list_of_intv=None):
        if list_of_intv is None:
            list_of_intv = self.admissible_gammas
        print(string)
        print('-------------- report ------------------')
        print('total width = {}'.format(self.get_total_width()))
        print('num intervals= {}'.format(len(list_of_intv)))
        for tn in list_of_intv:
            print(tn['interval'])
            print('samples {}'.format(tn['samples'].size))
            print('first & last samples {0[0]}, {0[1]}'.format((tn['samples'][0], tn['samples'][-1])))
            print('nan values {}'.format(sum(np.isnan(tn['samples']))))
        print('----------------------------------------')

    def gen_init_gammas(self, interval):
        return [{'interval': interval, 'samples': self.gen_sample_gammas(interval)}]

    def refine_admissible_gammas(self):
        # if self.verbose:
        #     self.report('entering refine_admissible_gammas()')
        copied_gammas = copy.deepcopy(self.admissible_gammas)

        for interval in copied_gammas:
            model_decision = self.decide(interval['samples'])
            # boolean array for nonzero decisions
            nonzero_decisions = ~np.ma.getmask(np.ma.masked_equal(model_decision, 0))
            # boolean array for incompatible decisions
            incompatible_decisions = model_decision != self.decision
            # set incompatible decisions to NaN
            interval['samples'][np.logical_and(nonzero_decisions, incompatible_decisions)] = np.nan
        all_gammas = np.concatenate(tuple(x['samples'] for x in copied_gammas), axis=0)
        tot_samples = all_gammas.size
        if self.verbose:
            self.report('work on copied gammas done', list_of_intv=copied_gammas)
        num_nan_samples = np.count_nonzero(np.isnan(all_gammas))
        if 0 < num_nan_samples < tot_samples:
            self.reconstruct_admissible_gammas(copied_gammas)
            self.total_width = self.get_total_width()
        elif num_nan_samples == tot_samples and self.verbose:
            print('depletion avoided')
        elif num_nan_samples == 0 and self.verbose:
            print('all samples are valid')
        self.refined = True

    def get_total_width(self):
        return sum([x['interval'][1] - x['interval'][0] for x in self.admissible_gammas])

    def reconstruct_admissible_gammas(self, intervals):
        interval_list = []
        for intvl in intervals:
            interval_list += self.reconstruct_interval(intvl)
        if not interval_list:  # if list is empty
            print('PROBLEM: INTERVAL RECONSTRUCTED IS EMPTY')
        self.admissible_gammas = [{'interval': x, 'samples': self.gen_sample_gammas(x)} for x in interval_list]
        if self.verbose:
            self.report('reconstruction of gammas over')

    def reconstruct_interval(self, interval):
        """
        splits an interval into sub-intervals
        :param interval: dict with: dict['interval'] = (lower_bound, upper_bound) and dict['samples'] = numpy.array
        :return: list of sub-intervals, possibly empty
        """
        interval_list = []
        last_lower_bound, last_upper_bound = interval['interval']
        old_samples = interval['samples']
        lower_bound = None  # new lower bound
        last_up = 0  # new upper bound

        for indx, g in enumerate(old_samples):
            nan = np.isnan(g)
            if (lower_bound is None) and nan:
                continue
            elif g <= last_up:
                continue
            elif lower_bound is None:
                nlb = g - self.tolerance_gamma
                if nlb >= max(last_lower_bound, last_up):
                    lower_bound = nlb
                else:
                    lower_bound = max(last_lower_bound, last_up)
                # for debug purposes
                if lower_bound < last_up:
                    print('WARNING: overlapping consecutive intervals')
            elif nan:
                gm1 = old_samples[max(indx - 1, 0)]
                last_up = gm1 + self.tolerance_gamma
                # following if just for debugging
                if last_up > lower_bound:
                    interval_list += [(lower_bound, last_up)]
                else:
                    print('WARNING: negative length interval!!!')
                lower_bound = None
            # list of samples ended without nan values, so inherit upper bound from prev. setting
            elif g == old_samples[-1]:
                last_up = last_upper_bound
                interval_list += [(lower_bound, last_up)]
        return interval_list

    def gen_sample_gammas(self, interval, max_range=50, max_samples=1000, min_distance=0.001):
        """
        samples uniformly in the interval
        :param interval: tuple with left value smaller than right value
        :param max_range: maximum range in which max_samples should be used
        :param max_samples: maximum number of samples to use if interval[1]-interval[0]<max_range
        :param min_distance: minimum distance between two consecutive samples, to avoid too many samples
        :return: numpy array of gamma samples + warning message if max_samples is too small for tolerance
        """
        # todo: return error if interval[1]-interval[0] <= 0
        curr_range = interval[1] - interval[0]
        if curr_range <= max_range:
            samples, step = np.linspace(interval[0], interval[1], max_samples, endpoint=False, retstep=True)
            if step > self.tolerance_gamma and self.verbose:
                print('WARNING: step = {0} while tolerance = {1}'.format(step, self.tolerance_gamma))
            elif step < min_distance:
                samples = np.arange(start=interval[0], stop=interval[1], step=min_distance)
        else:
            samples = np.arange(start=interval[0], stop=interval[1], step=self.tolerance_gamma)
            if self.verbose:
                fmt_string = 'WARNING TRIAL {}: sample size = {} while max_samples = {}'
                print(fmt_string.format(self.number, samples.size, max_samples))
        return samples

    def decide(self, gammas_array, init_cond=0):
        """
        computes decision using the linear model
        :param gammas_array: numpy n-by-1 array of linear discounting terms
        :param init_cond: initial value of accumulation variable
        :return: -1, 0 or 1 for left, undecided and right
        """
        y = init_cond
        gammas = gammas_array.reshape((-1, 1))

        # right train
        right_train = self.stimulus[1].reshape((1, -1))
        y += np.exp(gammas @ right_train).sum(axis=1)

        # left train
        left_train = self.stimulus[0].reshape((1, -1))
        y -= np.exp(gammas @ left_train).sum(axis=1)

        return np.sign(y)

    def stopping_criterion(self, stopping_width):
        """
        :param stopping_width: desired precision on inferred gamma. Must be greater than self.tolerance_gamma
        :return: (decision to stop, explanatory message)
        """
        if self.lost_gamma():
            msg = 'PROBLEM: true gamma has been lost'
            stop = True
        elif self.total_width < stopping_width and self.num_intervals == 1:
            msg = 'stopping criterion has been met'
            stop = True
        else:
            stop = False
            msg = ''

        return stop, msg

    def lost_gamma(self):
        """
        :return: True if true_gamma not in self.admissible_gammas anymore
        """
        found = False
        for intvl in self.admissible_gammas:
            if intvl['interval'][0] <= self.true_gamma <= intvl['interval'][1]:
                found = True
        return not found


"""
----------------------------GENERAL PURPOSE FUNCTIONS
"""


def get_lambda_high(lambda_low, s):
    """
    returns the high click rate, given the low click rate and S as input.
    :param lambda_low: low click rate
    :param s: S=(lambda_high-lambda_low)/sqrt(lambda_low+lambda_high)
    :return: value of lambda_high that fits with S and lambda_low
    """
    return (2 * lambda_low + s ** 2 + s * np.sqrt(s ** 2 + 8 * lambda_low)) / 2


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


"""
----------------------------FUNCTIONS FOR PARAMETER INFERENCE
"""


def run(num_trials, click_rates, true_gamma, interrogation_time, hazard_rate, init_range,
        stim_on_the_fly=True, verbose=False, independent_trials=False, global_gammas=None,
        report_full_list=False, report_widths=True, report_trials=None, file_name=None,
        run_number=None, swmr=False):
    """
    :param num_trials: number of trials in the run
    :param click_rates:
    :param true_gamma:
    :param interrogation_time:
    :param hazard_rate:
    :param init_range: init interval of admissible gammas
    :param stim_on_the_fly: reads from db if false
    :param verbose:
    :param independent_trials: if false, passes the admissible gammas recursively
    :param global_gammas:
    :param report_full_list:
    :param report_widths: for hists
    :param report_trials:
    :param file_name: file name for hfd5 db
    :param run_number: starts at 0
    :return:
    """
    # loop over trials to construct the trial-dependent list of admissible gammas
    if report_full_list:
        trial_list = []
    if report_widths:
        trial_widths = []
    # todo: throw error if num_trials < 1
    if num_trials < 1:
        print('num_trials should be >= 1')
        return
    if not stim_on_the_fly:
        # prepare DB
        f = h5py.File(file_name, 'r')
        group_name = build_group_name((click_rates[0], click_rates[1], hazard_rate, interrogation_time))
        data_trials = f[group_name + '/trials']
        data_info = f[group_name + '/trial_info']
        data_dec = f[group_name + '/decision_lin']
    for lll in range(num_trials):
        trial_nb = lll + 1
        if stim_on_the_fly:
            stim_train, _, _ = gen_stim(gen_cp(interrogation_time, hazard_rate), click_rates[0], click_rates[1],
                                     interrogation_time)
            # generate trial and decision
            if lll == 0 or independent_trials:
                trial = Trial(stim_train, true_gamma, trial_nb, verbose=verbose, init_gamma_range=init_range)
            else:
                trial = Trial(stim_train, true_gamma, trial_nb, init_gammas=global_gammas, verbose=verbose)
        else:
            # read trial from db
            row_nb = lll + num_trials * run_number
            stim_train = tuple(data_trials[row_nb, :2])
            decision_value = data_dec[row_nb]
            # generate trial and decision
            if lll == 0 or independent_trials:
                trial = Trial(stim_train, true_gamma, trial_nb, decision_datum=decision_value,
                              verbose=verbose, init_gamma_range=init_range)
            else:
                trial = Trial(stim_train, true_gamma, trial_nb, decision_datum=decision_value,
                              init_gammas=global_gammas, verbose=verbose)

        # test gamma samples and refine admissible interval
        trial.refine_admissible_gammas()

        # update global variables
        if report_full_list:
            trial_list += [trial]
        if report_widths:
            if trial_nb in report_trials:
                trial_widths += [(trial.total_width, trial.number)]
        global_gammas = copy.deepcopy(trial.admissible_gammas)

        # stopping criteria in addition to reaching num_trials in the upcoming for loop
        # exit for loop if stopping criterion met
        stop_loop, message = trial.stopping_criterion(1.01 * trial.tolerance_gamma)
        if stop_loop:
            print(message)
            break
    if not stim_on_the_fly:
        f.close()
    if report_full_list:
        return trial_list
    if report_widths:
        return trial_widths


def create_hfd5_data_structure(file, groupname, num_trials):
    """
    :param file: h5py.File
    :param groupname:
    :param num_trials: nb of trials
    :return: created group
    """
    group = file.create_group(groupname)
    dt = h5py.special_dtype(vlen=np.dtype('f'))
    group.create_dataset('trials', (num_trials, 3), maxshape=(100000, 10), dtype=dt)
    group.create_dataset('trial_info', (num_trials, 3), maxshape=(100000, 10), dtype='i')
    return group


def build_group_name(four_p):
    """
    :param four_p: (low click rate, high click rate, hazard rate, interrogation time)
    :return: string
    """
    return 'lr' + str(four_p[0]) + 'hr' + str(four_p[1]) + 'h' + str(four_p[2]) + 'T' + str(four_p[3])


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
        return 1.45333 + 0.670241 * snr + 0.34324 * (snr ** 2) - 0.00275835 * (snr ** 3)
    else:
        corr = {'gamma': np.array([2.0848, 2.5828, 3.3143, 4.2789, 5.4162, 6.7457, 8.1371, 9.7199,
                                   11.3937, 13.2381, 15.1327, 17.2771, 19.5909, 22.0435, 24.6947,
                                   27.7241, 30.5711, 33.5354, 36.7966, 40.3143]),
                'S/h': np.arange(0.5, 10.1, 0.5)}
        idx = np.where(corr['S/h'] == ratio_rates)[0][0]
        return corr['gamma'][idx]


def update_decision_data(file_name, model, group_name):
    """
    :param file_name: file name (string)
    :param model: string, either 'lin' or 'nonlin'
    :param group_name: group object from h5py module
    :return:
    """
    f = h5py.File(file_name, 'a')
    group = f[group_name]
    info_dset = group['trial_info']
    trials_dset = group['trials']
    num_trials = trials_dset.shape[0]
    dset_name = 'decision_' + model
    # check whether decision dataset exists
    if dset_name in group:
        dset = group[dset_name]
        dset_length = dset.shape[0]
        # resize for new data
        if dset_length < num_trials:
            dset.resize(num_trials, axis=0)
            row_indices = np.arange(dset_length, num_trials)
        else:
            return
    else:
        # create
        dset = group.create_dataset(dset_name, (num_trials, ), dtype='i', maxshape=(100000, ))
        row_indices = np.arange(num_trials)
    # store best gamma as attribute for future reference
    best_gamma = get_best_gamma(round(info_dset.attrs['S'], 4), 1)
    dset.attrs['best_gamma'] = best_gamma
    # populate dataset
    for row_idx in row_indices:
        # get trial object
        stim = tuple(trials_dset[row_idx, :2])
        trial = Trial(stim, best_gamma, row_idx + 1)
        dset[row_idx] = trial.decision
    f.flush()
    f.close()


def dump_info(four_parameters, s):
    print('S value: {}'.format(s))
    print('low click rate: {}'.format(four_parameters[0]))
    print('high click rate: {}'.format(four_parameters[1]))
    print('hazard rate: {}'.format(four_parameters[2]))
    print('interr. time: {}'.format(four_parameters[3]))


if __name__ == '__main__':
    # parameter vectors
    a_S = [0.5, 3, 8]  # S parameter
    # a_gamma = [2.0848, 6.7457, 27.7241]  # best gamma
    a_ll = [30, 15, 1]  # low click rate

    # scalar parameters
    int_time = 2
    hazard = 1

    S = a_S[1]
    lambda_low = a_ll[1]  # low click rate
    lambda_high = get_lambda_high(lambda_low, S)
    four_params = (lambda_low, lambda_high, hazard, int_time)
    true_g = get_best_gamma(S, hazard)
    start_time = time.time()

    filename = '/scratch/adrian/srvr_data_3.h5'
    # populate_hfd5_db(filename, four_params, 50000)
    update_decision_data(filename, 'lin', build_group_name(four_params))

    # create response datasets for best linear and nonlinear models
    lambda_low = a_ll[2]
    for S in a_S:
        four_params = (lambda_low, get_lambda_high(lambda_low, S), hazard, int_time)
        populate_hfd5_db(filename, four_params, 100000)
        update_decision_data(filename, 'lin', build_group_name(four_params))

  #  num_run = 500
   # number_of_trials = 40
   # report_nb = np.floor(np.linspace(1, number_of_trials, 4))
   # init_interval = (0, 40)
   # widths = [[] for _ in range(len(report_nb))]  # empty list of lists of total widths. One list per trial nb
   # for run_nb in range(num_run):
        # print('\n ///////////////////')
        # print('run {}'.format(run_nb + 1))
    #    sim_trials = run(number_of_trials, (lambda_low, lambda_high),
    #                     true_g, int_time, hazard, init_interval, run_number=run_nb,
    #                     verbose=False, stim_on_the_fly=False, file_name=filename,
    #                     report_full_list=False, report_widths=True, report_trials=report_nb)
    #    for sim_trial in sim_trials:
    #        tnb = sim_trial[1]
    #        for idxx, nb in enumerate(report_nb):
    #            if tnb == nb:
    #                widths[idxx].append(sim_trial[0])
    #dump_info(four_params, S)
    print("--- {} seconds ---".format(time.time() - start_time))
    #
    #for idx, ttt in enumerate(widths):
    #    plt.figure(figsize=(6, 2.8))
    #    # plt.subplot(len(report_nb), 1, idx + 1)
    #    plt.hist(ttt, bins='auto', density=True)
    #    plt.axvline(np.mean(ttt), color='r', linestyle='-', linewidth=2)
    #    plt.title('trial {}'.format(int(report_nb[idx])))
    #    plt.xlim(init_interval)
    #    plt.xlabel('total width')
    #    plt.ylabel('density')
    #    plt.tight_layout()
    # plt.show()
        # plt.savefig('/home/radillo/Pictures/simulations/PPPP_02{}.svg'.format(idx), bbox_inches='tight')
        # plt.savefig('/home/adrian/tosubmit_home/{}_{}.svg'.format(,idx), bbox_inches='tight')
    # if len(sys.argv) > 1:
    #     filename = 'report{}'.format(sys.argv[1])
    #     plt.savefig('/home/radillo/Pictures/simulations/{}.png'.format(filename), bbox_inches='tight')
