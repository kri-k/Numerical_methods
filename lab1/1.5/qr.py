import argparse
from libs.qr_eigenvalue_algorithm import *


def setup_argparse():
    parser = argparse.ArgumentParser(description=
                                     "Find eigenvalues (real and complex) of matrix")
    parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    parser.add_argument("-p", "--precision", action="store", default="0.01", metavar="value",
                        type=float,
                        help="precision of computation")
    parser.add_argument(dest="matrix_file", action="store", help="input matrix")
    return parser


if __name__ == "__main__":
    parser = setup_argparse()
    args = parser.parse_args()
    try:
        m = Matrix.parse(args.matrix_file)
    except FileNotFoundError as e:
        print(e.filename, e.strerror)
        sys.exit(1)
    qr_algo = QREigenvalueAlgorithm(m)
    qr_algo.log(args.verbose)
    values = qr_algo.solve(args.precision)
    for i, val in enumerate(values):
        if type(val) is complex and val.imag == 0:
            print("l{} = {}".format(i, val.real))
        else:
            print("l{} = {}".format(i, val))

