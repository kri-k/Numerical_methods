import argparse
from libs.tridiagonal_matrix_algorithm import *


def setup_argparse():
    parser = argparse.ArgumentParser(description=
                                     "Solve a system of linear equations with tridiagonal matrix algorithm (TDMA)")
    parser.add_argument("-v", "--verbose", action="store_true", help="increase output verbosity")
    parser.add_argument(dest="matrix_file", action="store", help="left part of linear system")
    parser.add_argument(dest="vector_file", action="store", help="right part of linear system")
    return parser


if __name__ == "__main__":
    parser = setup_argparse()
    args = parser.parse_args()
    mtrx = Matrix.parse(args.matrix_file)
    vec = Matrix.parse(args.vector_file)
    tdma = TDMA(mtrx, vec)
    tdma.log(args.verbose)
    print(tdma.solve())

