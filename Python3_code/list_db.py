import sys
import h5py
def vv(n):
    print(n)
if len(sys.argv) > 1:
    with h5py.File(sys.argv[1], 'r') as f:
        f.visit(vv)
