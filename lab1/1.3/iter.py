from libs.iterative_method import *
import argparse


def setup_argparse():
    parser = argparse.ArgumentParser(description=
                                     "Solve a system of linear equations with iterative method")
    parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    parser.add_argument("-p", "--precision", action="store", default="0.01", metavar="value",
                        type=float,
                        help="precision of computation")
    parser.add_argument(dest="matrix_file", action="store", help="left part of linear system")
    parser.add_argument(dest="vector_file", action="store", help="right part of linear system")
    return parser


if __name__ == "__main__":
    parser = setup_argparse()
    args = parser.parse_args()
    try:
        m = Matrix.parse(args.matrix_file)
        v = Matrix.parse(args.vector_file)
    except FileNotFoundError as e:
        print(e.filename, e.strerror)
        sys.exit(1)
    itr_m = IterativeMethod(m, v)
    itr_m.log(args.verbose)
    print(itr_m.solve(args.precision))
