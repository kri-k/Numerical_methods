import sys
import argparse
from libs.exponential_spline import ExpSpline


def setup_argparse():
    parser = argparse.ArgumentParser(description=
                                     'Find and plot the interpolating cubic spline')
    # parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
    parser.add_argument('-p', '--plot', action='store_true', help='plot spline using matplotlib')
    parser.add_argument(dest='points_file', action='store', help='file with points to interpolate')
    parser.add_argument("-s", "--step", action="store", default="100", metavar="value",
                        type=int,
                        help="number of steps to plot")
    return parser


if __name__ == "__main__":
    args = setup_argparse().parse_args()
    try:
        file = open(args.points_file)
    except FileNotFoundError as e:
        print(e.filename, e.strerror)
        sys.exit(1)
    px = []
    py = []
    px_find = []
    py_find = []
    for line in file:
        l = line.split()
        if len(l) == 0 or l[0] == '#':
            continue
        if l[0] == '?':
            px_find.append(float(l[1]))
        else:
            px.append(float(l[0]))
            py.append(float(l[1]))
    file.close()
    spl = ExpSpline(*zip(px, py))
    py_find = [spl(x) for x in px_find]
    for i, x in enumerate(px_find):
        print(x, py_find[i])
    if args.plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError as e:
            print("Can't import matplotlib :'''(")
            print(e)
            sys.exit(1)
        step_num = args.step
        step = abs(px[-1] - px[0]) / step_num
        x = [px[0]]
        y = [spl(x[0])]
        for i in range(1, step_num):
            x.append(x[i-1] + step)
            y.append(spl(x[-1]))
        x.append(px[-1])
        y.append(spl(x[-1]))
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Interpolating exponential spline for "{}"'.format(args.points_file))
        plt.grid(True)
        plt.plot(px, py, 'bo', x, y, '')
        plt.plot(px_find, py_find, 'ro')
        plt.show()
