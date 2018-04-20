"""
The aim of this script is to infer the discounting rate used to produce simulated data.
The script contains two types of functions.
    The ones used to simulate the data.
    The ones used to infer the parameter.

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
font = {'family': 'DejaVu Sans',
        'weight': 'bold',
        'size': 22}

matplotlib.rc('font', **font)
# from fractions import Fraction
# from math import gcd

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


def decide(click_trains, discounting_rate, init_cond=0):
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


if __name__ == '__main__':
    # todo: better deal with sample depletion
    # test code for single trial
    a_S = [.5, 3, 8]
    a_gamma = [2.0848, 6.7457, 27.7241]
    T = 2
    h = 1
    a_ll = [30, 15, 1]  # low click rate
    for jjj in range(len(a_ll)):
        ll = a_ll[jjj]
        for kkk in range(len(a_S)):
            S = a_S[kkk]
            gamma = a_gamma[kkk]
            num_gammas = 1000
            gamma_range = np.linspace(0, 50, num_gammas)
            num_trials = 500
            # all_gammas = np.zeros((num_gammas, num_trials))
            all_gammas = gamma_range.copy()
            for jj in range(num_trials - 1):
                all_gammas = np.c_[all_gammas, gamma_range]
            for trial in range(num_trials):
                gammas_backup = gamma_range.copy()
                if np.sum(np.isnan(gamma_range)) == num_gammas - 1:
                    break
                stim_train, last_envt_state = gen_stim(gen_cp(T, h), ll, get_lambda_high(ll, S), T)
                d = decide(stim_train, gamma)
                # range_gammas = cov_trains2poly(np.array([0.2857142857142857, 0.3333333333]), np.array([99/100])))
                for ii in range(gamma_range.size):
                    if not np.isnan(gamma_range[ii]):
                        model_dec = decide(stim_train, gamma_range[ii])
                        if model_dec == 0:
                            print('decision=0 for gamma = ' + str(gamma_range[ii]))
                            print('%i left clicks and %i right clicks' % (len(stim_train[0]),
                                                                          len(stim_train[1])))
                        if model_dec != d:
                            gamma_range[ii] = np.nan
                            all_gammas[ii, trial:] = np.nan
                if all(np.isnan(gamma_range)):  # to cope with sample depletion
                    gamma_range = gammas_backup.copy()
                    break

            print('%i valid values after %i trials' % (gamma_range.size - np.sum(np.isnan(gamma_range)), trial+1))
            print('true gamma = %f' % gamma)
            print('admissible gammas:')
            for ii in range(gamma_range.size):
                if not np.isnan(gamma_range[ii]):
                    print(gamma_range[ii])
            idx = 3 * jjj + kkk + 1
            plt.subplot(3, 3, idx)
            plt.plot(all_gammas[:, :trial].transpose(), 'b-')
            plt.ylabel('admissible gamma')
            plt.xlabel('Trial')
            title_string = 'S = %f; h_low = %i' % (S, ll)
            plt.title(title_string)

    plt.show()
