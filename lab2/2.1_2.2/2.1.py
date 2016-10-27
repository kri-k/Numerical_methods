from math import *
from libs.newtons_method import *
from libs.fixed_point_iteration import *


if __name__ == "__main__":
    print("=======================")
    print("Newtons Method:")
    print("f(x) = ln(1 + x) - 2*x^2 + 1 = 0")
    f = MathFunction(lambda x: log1p(x) - 2 * x**2 + 1,
                     lambda x: 1 / (x + 1) - 4 * x)
    solution1 = NewtonsMethod(f).solve(1e-10,
                                       [(0.4, 1.)],
                                       log=True)
    print("Solution =", solution1[0])
    print("f({}) = {}".format(solution1[0], f(solution1[0])))

    print()

    print("=======================")
    print("Fixed Point Iteration Method:")
    print("g(x) = x - lam*f(x) = x")
    solution2 = FixedPointIteration(f).solve(1e-10,
                                             [(0.4, 1.)],
                                             log=True)
    print("Solution =", solution2[0])
    print("f({}) = {}".format(solution2[0], f(solution2[0])))
