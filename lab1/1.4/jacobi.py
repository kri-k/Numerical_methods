import argparse
from libs.jacobi_eigenvalue_algorithm import *


def setup_argparse():
    parser = argparse.ArgumentParser(description=
                                     "Find eigenvalues and eigenvectors of matrix")
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
    jacobi_algo = JacobiEigenvalueAlgorithm(m)
    jacobi_algo.log(args.verbose)
    a, u = jacobi_algo.solve(args.precision)
    values = JacobiEigenvalueAlgorithm.split_to_values(a)
    vectors = JacobiEigenvalueAlgorithm.split_to_vectors(u)
    for i, val in enumerate(values):
        print("l{} = {}".format(i, val))
    print()
    for i, vec in enumerate(vectors):
        print("x{} = ".format(i))
        print(vec)

    if args.verbose:
        print("=========== Control check ===========")
        print("Check Ax = lx\n")
        for i in range(len(values)):
            print("Check value and vector #", i+1, sep='')
            print("Ax = ")
            print(m * vectors[i])
            print("lx = ")
            tmp_v = vectors[i].copy()
            tmp_v *= values[i]
            print(tmp_v)

        print("Checking orthogonality:\n")
        for i in range(len(vectors) - 1):
            print("(x{}, x{}) = {}".format(i+1, i+2, vectors[i].scalar_product(vectors[i+1])))
