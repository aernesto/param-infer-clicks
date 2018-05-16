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
from official_fcns import *

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
            curr_sq_err, curr_width = get_block_width(reference_dec, decision_data, all_sample_values, sp_tol, true_param)
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

    # compute squared error
    squared_error

    # compute total width

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
              'block_number': 100,
              'trial_number': 50,
              'model_to_fit': 'lin',
              'reference_model': 'lin'}
    params['high_rate'] = get_lambda_high(params['low_rate'], S)
    if params['S'] in np.arange(0.5, 10.1, 0.5) and params['hazard_rate'] == 1:
        pol = False
    else:
        pol = True
    params['best_gamma'] = get_best_gamma(params['S'], params['hazard_rate'], polyfit=pol)
    params['group_name'] = build_group_name((params['low_rate'],
                                             params['high_rate'],
                                             params['hazard_rate'],
                                             params['T']))
    MSE, AvgWidth = deter_fit(params)
