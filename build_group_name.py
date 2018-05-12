import numpy as np
import sys
from official_fcns import build_group_name

if __name__ == '__main__':
    if len(sys.argv) == 5:
        lr = float(sys.argv[1])
        hr = float(sys.argv[2])
        h = float(sys.argv[3])
        T = float(sys.argv[4])
        print(build_group_name((lr, hr, h, T)))
    else:
        raise ValueError('wrong nb of command line args')
