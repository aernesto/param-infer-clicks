import numpy as np
import sys

if __name__ == '__main__':
    if len(sys.argv) > 2:
        lambda_low = float(sys.argv[1])
        lambda_high = float(sys.argv[2])
        print(np.log(lambda_high/lambda_low))
