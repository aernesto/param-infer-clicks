"""
The aim of this script is to infer the discounting rate used to produce simulated data.
The script contains two types of functions.
    The ones used to simulate the data.
    The ones used to infer the parameter.

"""
import numpy as np
from fractions import Fraction
from math import gcd

"""
----------------------------GENERAL PURPOSE FUNCTIONS
"""


def lcm(numbers):
    """computes the least common multiple of a list of positive integers"""
    # check that all numbers are integers
    if all([isinstance(x, int) and x > 0 for x in numbers]):
        m = 1
        # compute lcm iteratively
        for n in numbers:
            m = int(n * m / gcd(n, m))
        return m
    else:
        print('some numbers are not positive integers')
        return None


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


"""
----------------------------FUNCTIONS USED TO INFER PARAMETERS
"""


# def sign_gamma(trains, gamma_range=(0, 80)):
#     """old function used to compute the decision at several gamma values"""
#     step = 0.01
#     gammas = np.r_[gamma_range[0]:gamma_range[1]+step:step]
#     init_sign = decide(trains, gammas[0])
#     last_sign = init_sign
#     switches = []
#     if gammas.size > 0:
#         for g in gammas:
#             new_sign = decide(trains, g)
#             if last_sign != new_sign:
#                 switches += [g]
#                 last_sign = new_sign
#         end_sign = new_sign
#     else:
#         print('Error: gamma_range should not be empty')
#         return None
#     return init_sign, end_sign, switches


def get_range_acceptable_gammas(trains, dec):
    """
    given the data from a single trial, returns the range of acceptable gammas
    the idea is to convert
    """
    left_clicks = [Fraction.from_float(x).limit_denominator(20) for x in trains[0]]
    right_clicks = [Fraction.from_float(x).limit_denominator(20) for x in trains[1]]

    denominators = [x.denominator for x in left_clicks] + [x.denominator for x in right_clicks]

    clicks_lcm = lcm(denominators)

    numerators = [(x.numerator * clicks_lcm / x.denominator, -1) for x in left_clicks] + \
                 [(x.numerator * clicks_lcm / x.denominator, 1) for x in right_clicks]

    numerators.sort(key=lambda tup: tup[0])  # sorts in place according to decreasing numerator value

    num_array = np.array(numerators)  # col 0 = nums, col 1 = idx
    if num_array.size == 0:
        powers = np.array([])
    else:
        max_power = int(num_array[-1][0])
        print(max_power)
        powers = np.zeros(max_power + 1)
        for i in range(len(num_array)):
            powers[int(num_array[i, 0])] = int(num_array[i, 1])
    powers = np.flip(powers, 0)
    return powers


if __name__ == '__main__':
    # test code for single trial
    S = .5
    gamma = 2.0848
    T = 2
    h = 1
    ll = 30  # low click rate
    stim_train, last_envt_state = gen_stim(gen_cp(T, h), ll, get_lambda_high(ll, S), T)
    d = decide(stim_train, gamma)
    range_gammas = get_range_acceptable_gammas(stim_train, d)
    if range_gammas is None:
        print('None')
    else:
        print(range_gammas)
