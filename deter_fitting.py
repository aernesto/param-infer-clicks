"""
The aim of this script is to infer the discounting & hazard rate used to produce simulated data.
"""
import numpy as np
import matplotlib
# matplotlib.use('Agg')  # required on server to forbid X-windows usage
import pickle
import matplotlib.pyplot as plt
import copy
import time
import sys
import h5py
from sympy import *

font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 15}

matplotlib.rc('font', **font)
# from fractions import Fraction
# from math import gcd


class Trial:  # todo: make sure trial_duration is passed correctly in rest of script
    def __init__(self, trial_number, init_gammas, init_h, tolerance, true_gamma, true_h):
        """
        :param true_gamma: true gamma with which decision data should be computed
        :param trial_number: number within external sequence
        :param init_gammas: initial list of gamma samples to try
        :param init_h: initial list of h samples to try
        :param tolerance: width to add to each side of bordering samples
        """
        self.number = trial_number
        self.true_gamma = true_gamma
        self.true_h = true_h
        self.tolerance = tolerance
        if (decision_datum is None) and linear:
            self.decision = self.decide(np.array([self.true_gamma]))  # -1 for left, 1 for right, 0 for undecided
        else:
            self.decision = decision_datum
        if nonlinear:
            self.nonlin_dec = self.nonlin_decide(hazard_rate)
        else:
            self.nonlin_dec = None
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

    def nonlin_decide(self, hh, init_cond=0):
        """
        :param hh: hazard rate to use for the decision
        :param init_cond: initial condition for ODE
        :return:
        """
        y = init_cond
        t = 0
        right_clicks_left = self.stimulus[1].size
        left_clicks_left = self.stimulus[0].size
        clicks_left = left_clicks_left + right_clicks_left

        right_clicks = list(self.stimulus[1])
        left_clicks = list(self.stimulus[0])
        while clicks_left > 0:
            if right_clicks_left and left_clicks_left:
                if right_clicks[0] < left_clicks[0]:
                    nxt_click = right_clicks[0]
                    right_clicks = right_clicks[1:]
                    right_clicks_left -= 1
                    dwell = nxt_click - t
                    y = end_point_nonlin(y, dwell, hh)
                    y += self.kappa
                    t = nxt_click
                elif right_clicks[0] > left_clicks[0]:
                    nxt_click = left_clicks[0]
                    left_clicks = left_clicks[1:]
                    left_clicks_left -= 1
                    dwell = nxt_click - t
                    y = end_point_nonlin(y, dwell, hh)
                    y -= self.kappa
                    t = nxt_click
                else:
                    right_clicks_left -= 1
                    left_clicks_left -= 1
                    right_clicks = right_clicks[1:]
                    left_clicks = left_clicks[1:]
            elif right_clicks_left:
                nxt_click = right_clicks[0]
                right_clicks = right_clicks[1:]
                right_clicks_left -= 1
                dwell = nxt_click - t
                y = end_point_nonlin(y, dwell, hh)
                y += self.kappa
                t = nxt_click
            elif left_clicks_left:
                nxt_click = left_clicks[0]
                left_clicks = left_clicks[1:]
                left_clicks_left -= 1
                dwell = nxt_click - t
                y = end_point_nonlin(y, dwell, hh)
                y -= self.kappa
                t = nxt_click
            clicks_left = left_clicks_left + right_clicks_left
        dwell = self.trial_duration - t
        y = end_point_nonlin(y, dwell, hh)
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


def end_point_nonlin(init_cond, end_time, hr):
    if init_cond == 0:
        return 0
    else:
        return float(2*mpmath.acoth(np.exp(2*hr*end_time)*mpmath.coth(init_cond/2)))


def get_lambda_high(lamb_low, s):
    """
    returns the high click rate, given the low click rate and S as input.
    :param lamb_low: low click rate
    :param s: S=(lambda_high-lambda_low)/sqrt(lambda_low+lambda_high)
    :return: value of lambda_high that fits with S and lambda_low
    """
    return (2 * lamb_low + s ** 2 + s * np.sqrt(s ** 2 + 8 * lamb_low)) / 2


def deter_fit(p):
    """
    :param p: a dict of parameters with the following keys:
        ['S','low_rate','high_rate','hazard_rate','T','best_gamma',
        'filename','samples_params','ntrials','group_name']
    :return: dict containing four lists of admissible intervals
    """
    # prepare DB
    db_file = h5py.File(p['filename'], 'r')
    group_name = p['group_name']
    dset_trials = db_file[group_name + '/trials']
    dset_info = db_file[group_name + '/trial_info']
    dset_lin_dec = db_file[group_name + '/decision_lin']
    dset_nonlin_dec = db_file[group_name + '/decision_nonlin']
    T = p['T']
    fit_results_dict = {'lin2lin': [],
                        'lin2nonlin': [],
                        'nonlin2lin': [],
                        'nonlin2nonlin': []}
    all_sample_values = build_sample_vec(p['samples_params'])
    for model_condition in fit_results_dict.keys():

        fit_results_dict[model_condition] = get_intervals(model_condition)

    db_file.close()
    return fit_results_dict


def get_intervals(ref_dec, synthetic_dec, init_sample_values):
    """
    :param ref_dec: value of reference decision
    :param synthetic_dec: vector of decision per sample value
    :param init_sample_values: initial numpy array of admissible samples
    :return: list of admissible intervals
    """





def build_sample_vec(samples_params_dict):
    start = samples_params_dict['start']
    end = samples_params_dict['end']
    nb = samples_params_dict['number']
    return np.linspace(start, end, nb)

def build_group_name(four_p):
    """
    :param four_p: (low click rate, high click rate, hazard rate, interrogation time)
    :return: string
    """
    return 'lr' + str(four_p[0]) + 'hr' + str(four_p[1]) + 'h' + str(four_p[2]) + 'T' + str(four_p[3])


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
        iddx = np.where(corr['S/sqrt(h)'] == skellam/np.sqrt(h))[0][0]
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
              'low_rate': 1,
              'hazard_rate': 1,
              'T': 2,
              'filename': 'data/',
              'samples_params': {'start': 0, 'end': 40, 'number': 10000},
              'ntrials': 100000}
    params['high_rate'] = get_lambda_high(params['low_rate'], S)
    params['best_gamma'] = get_best_gamma(params['S'], params['hazard_rate'])
    params['group_name'] = build_group_name((params['low_rate'],
                                             params['high_rate'],
                                             params['hazard_rate'],
                                             params['T']))
    admissible_intervals = deter_fit(params)
