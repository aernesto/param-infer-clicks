"""
The aim of this script is to infer the discounting rate used to produce simulated data.
The script contains two types of functions.
    The ones used to simulate the data.
    The ones used to infer the parameter.

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import copy
font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 22}

matplotlib.rc('font', **font)
# from fractions import Fraction
# from math import gcd

"""
----------------------------CLASS DEFINITIONS
"""

# class Gammas(list):
    # def __init__(self, num_trials=1):
    #     self.num_trials = num_trials


class Trial:
    def __init__(self, stimulus, gamma, init_gammas=None, tolerance=.05):
        """
        :param stimulus: 2-tuple of click trains (left, right)
        :param gamma: true gamma with which decision data should be computed
        :param init_gammas: initial intervals of admissible gammas (list of dicts)
        :param tolerance: max distance allowed between two consecutive samples
        """
        self.stimulus = stimulus
        self.true_gamma = gamma
        self.tolerance_gamma = tolerance
        self.decision = self.decide(self.stimulus, self.true_gamma)  # -1 for left, 1 for right, 0 for undecided
        if init_gammas is None:
            self.admissible_gammas = self.gen_init_gammas()
        else:
            self.admissible_gammas = init_gammas
        self.refined = False  # True if self.refine_admissible_gammas() has been called
        self.total_width = self.get_total_width()
        self.num_intervals = len(self.admissible_gammas)

    def gen_init_gammas(self, interval=(0, 50)):
        return [{'interval': interval, 'samples': self.gen_sample_gammas(interval)}]

    def refine_admissible_gammas(self):
        copied_gammas = copy.deepcopy(self.admissible_gammas)
        for interval in copied_gammas:
            for idx, gamma in enumerate(interval['samples']):
                model_decision = self.decide(self.stimulus, gamma)
                if (model_decision != self.decision) and not model_decision:
                    interval['samples'][idx] = np.nan
        all_gammas = np.concatenate(tuple(x['samples'] for x in copied_gammas), axis=0)
        # print(type(all_gammas))
        # print('size all gammas = %i' % all_gammas.size)
        if not all(np.isnan(all_gammas)):
            self.reconstruct_admissible_gammas(all_gammas)
            # if all_gammas.size == 1:
            #     self.reconstruct_admissible_gammas([all_gammas])
            # else:
            #     self.reconstruct_admissible_gammas(all_gammas)
        else:
            print('depletion avoided')
        self.total_width = self.get_total_width()
        self.refined = True

    def get_total_width(self):
        return sum([x['interval'][1]-x['interval'][0] for x in self.admissible_gammas])

    def reconstruct_admissible_gammas(self, allg):
        # print('entering fcn')
        # print(type(allg))
        # print('first 5 elements of allg are')
        # print(allg[:5])
        # extract intervals first
        interval_list = []
        lower_bound = None
        last_up = 0

        # find lower_bound
        # find upper_bound
        # add interval to list
        it = np.nditer(allg, flags=['f_index'])
        while not it.finished:
            idx = it.index
            g = it[0]
            # print(it.index)
            # print(it[0])
            # print('-----')
            # print('idx %i' % idx)
            # print('gamma %f' % g)
            nan = np.isnan(g)
            it.iternext()
            if (lower_bound is None) and nan:
                continue
            elif lower_bound is None:
                nlb = g - self.tolerance_gamma
                if nlb > last_up:
                    lower_bound = nlb
                else:
                    continue
            elif nan:
                gm1 = allg[idx - 1]
                last_up = gm1 + self.tolerance_gamma
                interval_list += [(lower_bound, last_up)]
                lower_bound = None

        self.admissible_gammas = [{'interval': x, 'samples': self.gen_sample_gammas(x)} for x in interval_list]

    def gen_sample_gammas(self, interval, max_samples=1000):
        """
        samples uniformly in the interval
        :param interval: tuple with left value smaller than right value
        :param max_samples: maximum number of samples
        :return: numpy array of gamma samples + warning message if max_samples is too small for tolerance
        """
        samples = np.arange(start=interval[0], stop=interval[1], step=self.tolerance_gamma)
        if samples.size > max_samples:
            print('WARNING: one of tolerance or max_samples is too small. Resampling according to max_samples')
            samples, step = np.linspace(interval[0], interval[1], max_samples, endpoint=False, retstep=True)
            print('new between-samples step = %f' % step)
        return samples

    def decide(self, click_trains, discounting_rate, init_cond=0):
        """
        computes decision, given input trains
        :param click_trains: tuple of left and right click streams (numpy arrays)
        :param discounting_rate: linear discounting term
        :param init_cond: initial value of accumulation variable
        :return: -1, 0 or 1 for left, undecided and right
        """
        y = init_cond
        # right train
        for i in click_trains[1]:
            y += np.exp(discounting_rate * i)
        # left train
        for j in click_trains[0]:
            y -= np.exp(discounting_rate * j)
        return np.sign(y)

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

    return (np.array(left_stream), np.array(right_stream)), state[-1]


def stopping_criterion(width, ninter, stopping_width=0.01):
    """
    :param width: current precision on inferred gamma
    :param ninter: number of retained admissible gamma intervals
    :param stopping_width: desired precision on inferred gamma
    :return: True or False
    """
    return (width < stopping_width) and (ninter == 1)


if __name__ == '__main__':
    # todo: measure computation time
    # test code for single trial
    a_S = [.5]  # , 3, 8]
    a_gamma = [2.0848]  # , 6.7457, 27.7241]
    T = 2
    h = 1
    a_ll = [30, 15, 1]  # low click rate
    init_interval = (0, 50)  # initial interval of admissible gammas
    num_trials = 3
    for jjj in range(len(a_ll)):
        ll = a_ll[jjj]
        for kkk in range(len(a_S)):
            S = a_S[kkk]
            true_gamma = a_gamma[kkk]
            # num_gammas = 1000  # max sample size per admissible interval
            # all_gammas = np.zeros((num_gammas, num_trials))
            # all_gammas = gamma_range.copy()
            # for jj in range(num_trials - 1):
            #     all_gammas = np.c_[all_gammas, gamma_range]

            # loop over trials to construct the trial-dependent list of admissible gammas
            trial_list = []
            for lll in range(num_trials):
                stim_train, _ = gen_stim(gen_cp(T, h), ll, get_lambda_high(ll, S), T)

                # generate trial and decision
                if lll == 0:
                    trial = Trial(stim_train, true_gamma)
                else:
                    trial = Trial(stim_train, true_gamma, init_gammas=global_gammas)

                # test gamma samples and refine admissible interval
                trial.refine_admissible_gammas()

                # update global variables
                trial_list += [trial]
                global_gammas = copy.deepcopy(trial.admissible_gammas)

                # stopping criteria in addition to reaching num_trials in the upcoming for loop
                # exit for loop if stopping criterion met
                if stopping_criterion(trial.total_width, trial.num_intervals):
                    print('stopping criterion met')
                    break
            # gamma_samples = global_gammas[0]['samples']
            # print('%i valid values after %i trials' % (global_gammas[0]['samples'].size, lll+1))
            # print('true gamma = %f' % true_gamma)
            # print('admissible gammas:')
            # for ii in range(gamma_samples.size):
            #     if not np.isnan(gamma_samples[ii]):
            #         print(gamma_samples[ii])

            # idx = 3 * jjj + kkk + 1
            # plt.subplot(3, 3, idx)
            # # plt.plot(all_gammas[:, :trial].transpose(), 'b-')
            # plt.ylabel('admissible gamma')
            # plt.xlabel('Trial')
            # title_string = 'S = %f; h_low = %i' % (S, ll)
            # plt.title(title_string)

    # plt.show()
