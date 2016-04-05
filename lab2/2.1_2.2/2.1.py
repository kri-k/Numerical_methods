from math import *
from libs.newtons_method import *
from libs.fixed_point_iteration import *


if __name__ == "__main__":
    print("=======================")
    print("Newtons Method:")
    print("f(x) = ln(1 + x) - 2*x^2 + 1 = 0")
    f1 = MathFunction(lambda x: log1p(x) - 2 * x**2 + 1,
                      lambda x: 1 / (x + 1) - 4 * x)
    solution1 = NewtonsMethod(f1).solve(1e-10,
                                        (0.5,),
                                        log=True)
    print("Solution =", solution1[0])
    print("f({}) = {}".format(solution1[0], f1(solution1[0])))

    print()

    print("=======================")
    print("Fixed Point Iteration Method:")
    print("g(x) = sqrt((ln(1 + x) + 1) / 2) = x")
    f2 = MathFunction(lambda x: sqrt((log1p(x) + 1) / 2))
    solution2 = FixedPointIteration(f2).solve(1e-10,
                                              (0.5,),
                                              (0.45,),
                                              log=True)
    print("Solution =", solution2[0])
    print("f({}) = {}".format(solution2[0], f1(solution2[0])))
