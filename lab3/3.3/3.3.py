import sys
import argparse
from libs.ols import OLS


def setup_argparse():
    parser = argparse.ArgumentParser(description=
                                     'Find and plot the interpolating cubic spline')
    # parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
    parser.add_argument('-p', '--plot', action='store_true', help='plot spline using matplotlib')
    parser.add_argument('-e', '--error', action='store_true', help='compute error sum of squares')
    parser.add_argument(dest='points_file', action='store', help='file with points to approximate')
    parser.add_argument("-s", "--step", action="store", default="100", metavar="value",
                        type=int,
                        help="number of steps to plot")
    return parser


if __name__ == "__main__":
    parser = setup_argparse()
    args = parser.parse_args()
    try:
        file = open(args.points_file)
    except FileNotFoundError as e:
        print(e.filename, e.strerror)
        sys.exit(1)
    px = []
    py = []
    orders = []
    for line in file:
        l = line.split()
        if len(l) == 0 or l[0] == '#':
            continue
        if l[0].lower() == 'order':
            orders += map(int, l[1:])
        else:
            px.append(float(l[0]))
            py.append(float(l[1]))
    file.close()
    ols_list = [OLS(*zip(px, py), order=i) for i in orders]
    if args.error:
        print('Error sum of squares:')
        print('Order\tError')
        for i, x in enumerate(orders):
            print(x, ols_list[i].error_sum_squares(), sep='\t')
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
        y_orders_lists = [[f(x[-1])] for f in ols_list]
        for i in range(1, step_num):
            x.append(x[i-1] + step)
            for j, f in enumerate(ols_list):
                y_orders_lists[j].append(f(x[-1]))
        x.append(px[-1])
        for j, f in enumerate(ols_list):
            y_orders_lists[j].append(f(x[-1]))

        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Approximating with Ordinary Least Squares (OLS) for "{}"'.format(args.points_file))
        plt.grid(True)
        plt.plot(px, py, 'ro')
        for y in y_orders_lists:
            plt.plot(x, y)
        plt.show()
