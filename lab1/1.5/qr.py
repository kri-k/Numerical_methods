import argparse
from libs.qr_eigenvalue_algorithm import *


def setup_argparse():
    parser = argparse.ArgumentParser(description=
                                     "Find eigenvalues (real and complex) of matrix")
    parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    parser.add_argument("-s", "--step", action="store_true", help="show step by step solution")
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
    qr_algo.step_by_step(args.step)
    values = qr_algo.solve(args.precision)
    for i, val in enumerate(values):
        print("l{} = {}".format(i, val))
    print()
    val_vec = qr_algo.find_eigenvectors(values, 0.001, args.precision)
    for i, x in enumerate(val_vec):
        print("x{} =".format(i))
        print(x[1])
    for i, x in enumerate(val_vec):
        print("========================")
        a_x = m * x[1]
        l_x = x[1].copy()
        l_x *= x[0]
        print("A*x{} = ".format(i))
        print(a_x)
        print("l{0}*x{0} = ".format(i))
        print(l_x)
        print("delta = ")
        print(a_x - l_x)
