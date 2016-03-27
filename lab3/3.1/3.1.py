from math import *
from libs.lagrange_polynomial import LagrangePolynomial
from libs.newton_polynomial import NewtonPolynomial


def do_stuff(true_func,
             interpolated_func,
             points,
             x_check,
             description,
             func_name):
    print("==================")
    print(description)
    print(*points, sep=", ")
    print("{}({}) = {}".format(func_name, x_check, interpolated_func(x_check)))
    print("Error of interpolation =", abs(true_func(x_check) - interpolated_func(x_check)))


if __name__ == "__main__":
    f = lambda x: asin(x)
    x_a = [-0.4, -0.1, 0.2, 0.5]
    x_b = [-0.4, 0., 0.2, 0.5]
    points_a = [(x, f(x)) for x in x_a]
    points_b = [(x, f(x)) for x in x_b]
    lp_a = LagrangePolynomial(*points_a)
    lp_b = LagrangePolynomial(*points_b)
    np_a = NewtonPolynomial(*points_a)
    np_b = NewtonPolynomial(*points_b)
    x_check = 0.1
    print("f(x) = arcsin(x)")
    print("f({}) = {}".format(x_check, f(x_check)))
    do_stuff(f, lp_a, points_a, x_check, "Lagrange Polynomial lp_a(x) for points:", "lp_a")
    do_stuff(f, lp_b, points_b, x_check, "Lagrange Polynomial lp_b(x) for points:", "lp_b")
    do_stuff(f, np_a, points_a, x_check, "Newton Polynomial np_a(x) for points:", "np_a")
    do_stuff(f, np_b, points_b, x_check, "Newton Polynomial np_b(x) for points:", "np_b")

