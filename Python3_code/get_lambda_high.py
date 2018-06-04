import numpy as np
import sys
from official_fcns import get_lambda_high

if __name__ == '__main__':
    if len(sys.argv) > 2:
        lambda_low = float(sys.argv[1])
        S = float(sys.argv[2])
        print(get_lambda_high(lambda_low, S))
