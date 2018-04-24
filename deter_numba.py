"""
The aim of this script is to infer the discounting rate used to produce simulated data.
The script contains two types of functions.
    The ones used to simulate the data.
    The ones used to infer the parameter.

"""
import numpy as np
import matplotlib
# matplotlib.use('Agg')  # required on server to forbid X-windows usage
import matplotlib.pyplot as plt
import copy
import time
import sys
from numba import jit
from numba import jitclass          # import the decorator
from numba import int32, float32    # import the types

font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 16}

matplotlib.rc('font', **font)
# from fractions import Fraction
# from math import gcd

"""
----------------------------CLASS DEFINITIONS
"""


class Trial:
    def __init__(self, stimulus, gamma, trial_number, init_gammas=None, tolerance=.05, verbose=False):
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
        self.decision = self.decide(self.true_gamma)  # -1 for left, 1 for right, 0 for undecided
        if init_gammas is None:
            self.admissible_gammas = self.gen_init_gammas()
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

    def gen_init_gammas(self, interval=(0, 50)):
        return [{'interval': interval, 'samples': self.gen_sample_gammas(interval)}]

    def refine_admissible_gammas(self):
        # if self.verbose:
        #     self.report('entering refine_admissible_gammas()')
        copied_gammas = copy.deepcopy(self.admissible_gammas)

        for interval in copied_gammas:
            for indx, gamma in enumerate(interval['samples']):
                model_decision = self.decide(gamma)
                if (model_decision != self.decision) and model_decision:
                    interval['samples'][indx] = np.nan
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
        return sum([x['interval'][1]-x['interval'][0] for x in self.admissible_gammas])

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

    def decide(self, discounting_rate, init_cond=0):
        """
        computes decision, given discounting rate
        :param discounting_rate: linear discounting term
        :param init_cond: initial value of accumulation variable
        :return: -1, 0 or 1 for left, undecided and right
        """
        y = init_cond
        # right train
        for i in self.stimulus[1]:
            y += np.exp(discounting_rate * i)
        # left train
        for j in self.stimulus[0]:
            y -= np.exp(discounting_rate * j)
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


@jit
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


@jit
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


@jit
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

    return (np.array(left_stream), np.array(right_stream)), state[-1]


"""
----------------------------FUNCTIONS FOR PARAMETER INFERENCE
"""


@jit
def run(num_trials, click_rates, true_gamma, interrogation_time, hazard, stim_on_the_fly=True, verbose=False,
        independent_trials=False, global_gammas=None):
    # loop over trials to construct the trial-dependent list of admissible gammas
    trial_list = []
    # todo: throw error if num_trials < 1
    for lll in range(num_trials):
        trial_nb = lll + 1
        if stim_on_the_fly:
            stim_train, _ = gen_stim(gen_cp(interrogation_time, hazard), click_rates[0], click_rates[1],
                                     interrogation_time)

        # generate trial and decision
        if lll == 0 or independent_trials:
            trial = Trial(stim_train, true_gamma, trial_nb, verbose=verbose)
        else:
            trial = Trial(stim_train, true_gamma, trial_nb, init_gammas=global_gammas, verbose=verbose)

        # test gamma samples and refine admissible interval
        trial.refine_admissible_gammas()

        # update global variables
        trial_list += [trial]
        global_gammas = copy.deepcopy(trial.admissible_gammas)

        # stopping criteria in addition to reaching num_trials in the upcoming for loop
        # exit for loop if stopping criterion met
        stop_loop, message = trial.stopping_criterion(1.01 * trial.tolerance_gamma)
        if stop_loop:
            print(message)
            print("--- %s seconds ---" % (time.time() - start_time))
            break
        if lll == num_trials - 1:
            print('all trials used for refinement')
    return trial_list


if __name__ == '__main__':
    # test code for single trial
    a_S = [.5]  # 3, 8]
    a_gamma = [2.0848]  # 6.7457, 27.7241]  # best gamma
    T = 2
    h = 1
    a_ll = [30, 15, 1]  # low click rate
    init_interval = (0, 50)  # initial interval of admissible gammas
    number_of_trials = 6
    for jjj in [0]:  # range(len(a_ll)):
        ll = a_ll[jjj]
        for kkk in [0]:  # range(len(a_S)):
            start_time = time.time()

            S = a_S[kkk]
            true_g = a_gamma[kkk]
            lh = get_lambda_high(ll, S)
            num_run = 2
            report_nb = [1, 10, 20]
            widths = [[] for _ in range(len(report_nb))]  # empty list of lists of total widths. One list per trial nb
            for run_nb in range(num_run):
                print('\n ///////////////////')
                print('run {}'.format(run_nb+1))
                sim_trials = run(number_of_trials, (ll, lh), true_g, T, h, verbose=False)
                for tt in sim_trials:
                    tnb = tt.number
                    for idxx, nb in enumerate(report_nb):
                        if tnb == nb:
                            widths[idxx].append(tt.total_width)
            print("--- {} seconds ---".format(time.time() - start_time))
            for idx, ttt in enumerate(widths):
                plt.subplot(3, 1, idx + 1)
                plt.hist(ttt)
                plt.title('trial {}'.format(report_nb[idx]))

            # idx = 3 * jjj + kkk + 1
            # plt.subplot(3, 3, idx)
            # trial_num = 0
            # for t in sim_trials:
            #     trial_num += 1
            #     for g_intvl in t.admissible_gammas:
            #         plt.plot([trial_num, trial_num], [g_intvl['interval'][0], g_intvl['interval'][1]], 'b-')
            # plt.plot([1, len(sim_trials)], [true_g, true_g], 'r-')
            # plt.ylabel('admissible gamma')
            # plt.xlabel('Trial')
            # title_string = 'S = %f; h_low = %i' % (S, ll)
            # plt.title(title_string)

    # plt.show()
    # plt.savefig('/scratch/adrian/HISTS.png', bbox_inches='tight')
    if len(sys.argv) > 1:
        filename = 'report{}'.format(sys.argv[1])
        plt.savefig('/home/radillo/Pictures/simulations/{}.png'.format(filename), bbox_inches='tight')