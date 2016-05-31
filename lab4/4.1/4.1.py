import math
import matplotlib.pyplot as plt
from libs.сauchy_problem import CauchyProblem


def test():
    """
    var #8

    y'' - 4xy' + (4x^2 - 3)y - e^x^2 = 0,
    y(0) = 1,
    y'(0) = 0,
    x in [0, 1], h = 0.1

    Transform it to the system:
    y' = z = f1(x, y, z),
    z' = 4xz - (4x^2 - 3)y + e^x^2 = f2(x, y, z),
    y(0) = 1,
    z(0) = 0,
    x in [0, 1], h = 0.1

    Analytical solution is:
    y = (e^x + e^-x - 1)e^x^2

    """
    interval = [0, 1]
    step = 0.1

    def plot(sol, lbl):
        x = [x[0] for x in sol]
        y = [x[1][0] for x in sol]
        plt.plot(x, y, label=lbl)

    problem = CauchyProblem([1, 0],
                            lambda x, y, z: z,
                            lambda x, y, z: 4*x*z - (4*x**2 - 3)*y + math.e**x**2)
    methods = ('euler', 'runge_kutta', 'adams')

    solutions = [problem.solve(interval, step, m) for m in methods]
    rr_solutions = [CauchyProblem.runge_romberg(problem, interval, step, m) for m in methods]
    analytical_solution = [(x[0], [(math.e**x[0] + math.e**-x[0] - 1)*math.e**x[0]**2]) for x in solutions[0]]

    plot(analytical_solution, 'analytical solution')
    for i in range(len(methods)):
        plot(solutions[i], methods[i])
    plt.legend()

    for i in range(len(methods)):
        plt.figure()
        plot(analytical_solution, 'analytical solution')
        plot(solutions[i], methods[i])
        plot(rr_solutions[i], 'runge romberg + ' + methods[i])
        plt.legend()

    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

test()
