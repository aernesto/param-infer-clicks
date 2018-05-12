import numpy as np
import sys


def get_lambda_high(lamb_low, s):
    """
    returns the high click rate, given the low click rate and S as input.
    :param lamb_low: low click rate
    :param s: S=(lambda_high-lambda_low)/sqrt(lambda_low+lambda_high)
    :return: value of lambda_high that fits with S and lambda_low
    """
    return (2 * lamb_low + s ** 2 + s * np.sqrt(s ** 2 + 8 * lamb_low)) / 2


if __name__ == '__main__':
    if len(sys.argv) > 2:
        lambda_low = float(sys.argv[1])
        S = float(sys.argv[2])
        print(get_lambda_high(lambda_low, S))
