import argparse
from libs.gauss_method import *


def setup_argparse():
    parser = argparse.ArgumentParser(description=
                                     'Solve a system of linear equations with Gauss method and LUP decomposition')
    parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
    tpl = ('left_part_file', 'right_part_file')
    parser.add_argument('-s', '--solve', nargs=2, metavar=tpl, help='solve a system')
    parser.add_argument('-i', '--inverse', metavar='matrix_file', help='find the inverse matrix')
    parser.add_argument('-d', '--determinant',  metavar='matrix_file', help='find the determinant of a matrix')
    return parser


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Type '{} -h' for help".format(os.path.split(sys.argv[0])[-1]))
    parser = setup_argparse()
    args = parser.parse_args()
    if args.solve is not None:
        m = Matrix.parse(args.solve[0])
        v = Matrix.parse(args.solve[1])
        gm = GaussMethod(m, v)
        gm.log(args.verbose)
        print(gm.solve())
    if args.inverse is not None:
        m = Matrix.parse(args.inverse)
        print(GaussMethod.get_inverse_matrix(m, args.verbose))
    if args.determinant is not None:
        m = Matrix.parse(args.determinant)
        print(GaussMethod.get_determinant(m, args.verbose))
