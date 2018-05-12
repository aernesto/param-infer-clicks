import numpy as np
import sys

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
