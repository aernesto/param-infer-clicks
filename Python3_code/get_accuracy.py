import sys
from official_fcns import get_accuracy

if __name__ == '__main__':
    if len(sys.argv) > 3:
        filename = sys.argv[1]
        grpname = sys.argv[2]
        model = sys.argv[3]
        try:
            splecol = int(sys.argv[4])
        except IndexError:
            splecol = 0
        print(get_accuracy(filename, grpname, model, splecol))
