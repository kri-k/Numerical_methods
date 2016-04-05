from math import *
from libs.newtons_method import *
from libs.fixed_point_iteration import *


if __name__ == "__main__":
    print("=======================")
    print("Newtons Method:")
    print("f1(x1, x2) = x1^2 + x2^2 - 3^2 = 0")
    print("f2(x1, x2) = x1 - e^x2 + 3 = 0")
    f1 = MathFunction(lambda x1, x2: x1**2 + x2**2 - 3**2,
                      lambda x1, x2: 2*x1,
                      lambda x1, x2: 2*x2)
    f2 = MathFunction(lambda x1, x2: x1 - e**x2 + 3,
                      lambda x1, x2: 1,
                      lambda x1, x2: -e**x2)
    solution1 = SystemNewtonsMethod(f1, f2).solve(1e-10, [(2, 2)], log=True)
    print("x* =", solution1[0])
    print("f1(x*) = {}".format(f1(*solution1[0])))
    print("f2(x*) = {}".format(f2(*solution1[0])))

    print()

    print("=======================")
    print("Fixed Point Iteration Method:")
    print("g1(x1, x2) = sqrt(-x2^2 + 3^2) = x1")
    print("g2(x1, x2) = ln(x1 + 3) = x2")
    phi1 = MathFunction(lambda x1, x2: sqrt(-x2**2 + 3**2))
    phi2 = MathFunction(lambda x1, x2: log(x1 + 3))
    solution2 = SystemFixedPointIteration(phi1, phi2).solve(1e-10, [(2.5, 2.5)], [0.5], log=True)
    print("x* =", solution2[0])
    print("f1(x*) = {}".format(f1(*solution1[0])))
    print("f2(x*) = {}".format(f2(*solution1[0])))
