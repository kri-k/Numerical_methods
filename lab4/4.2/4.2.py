import matplotlib.pyplot as plt
from libs.boundary_value_problem import ShootingMethod, FiniteDifferenceMethod


def plot(sol, lbl):
    x = [x[0] for x in sol]
    y = [x[1][0] if type(x[1]) is list else x[1] for x in sol]
    plt.plot(x, y, label=lbl)


def test_1():
    """
    (x^2 + 1)y'' - 2y = 0,
    y'(0) = 0,
    y(2) - y'(2) = 1

    Transform to the system:
    y' = z,
    z' = 2y / (x^2 + 1),
    x in [0, 2],
    y(0) = n (parameter n varies),
    y'(0) = 0,
    condition of finishing calculations:
    y(2) - y'(2) - 1 < eps

    Analytical solution:
    y(x) = x**2 + 1
    """
    step = 0.2
    eps = 1e-5
    interval = [0, 2]
    pr = ShootingMethod(interval=interval,
                        initial_value_funcs=[lambda n: n, lambda n: 0],
                        functions=[lambda x, y, z: z, lambda x, y, z: 2*y / (x**2 + 1)],
                        right_boundary_func=lambda sol_b: sol_b[0] - sol_b[1] - 1)
    s = pr.solve(step, [5, 10], eps)
    real_s = [(x[0], x[0]**2 + 1) for x in s]
    rr_sol = pr.runge_romberg(step, [5, 10], eps)

    plt.title('step = %g' % step)
    plot(real_s, 'analytic solution')
    plot(s, 'shooting method')
    plot(rr_sol, 'runge-romberg')
    plt.legend()


def test_2():
    """
    (x^2 + 1)y'' - 2y = 0,
    y'(0) = 0,
    y(2) - y'(2) = 1

    (x^2 + 1)y'' - 2y = 0 transform to:
    y'' - 2 / (x^2 + 1) * y = 0
    =>
    p(x) = 0
    q(x) = 2 / (x^2 + 1)
    f(x) = 0

    Analytical solution:
    y(x) = x**2 + 1
    """
    def get_row_1(pr_obj, x, h):
        a, b, c, d = pr_obj.get_coef_tuple(x, h)
        return a + b, c, d

    def get_row_n(pr_obj, x, h):
        x -= h
        a, b, c, d = pr_obj.get_coef_tuple(x, h)
        return a, b - c / (h-1), d - c * h / (h - 1)

    p = lambda x: 0
    q = lambda x: -2 / (x**2 + 1)
    f = lambda x: 0

    step = 0.2
    interval = [0, 2]
    pr = FiniteDifferenceMethod(interval, p, q, f)
    s = pr.solve(step, get_row_1(pr, interval[0], step), get_row_n(pr, interval[1], step),
                 left_value=lambda y, h: y,
                 right_value=lambda y, h: (h - y) / (h - 1))
    real_s = [(x[0], x[0]**2 + 1) for x in s]
    rr_sol = pr.runge_romberg(step,
                              get_row_1,
                              get_row_n,
                              left_value=lambda y, h: y,
                              right_value=lambda y, h: (h - y) / (h - 1))
    plt.title('step = %g' % step)
    plot(real_s, 'analytic solution')
    plot(s, 'finite difference solution')
    plot(rr_sol, 'runge-romberg')
    plt.legend()

test_1()
plt.figure()
test_2()
plt.show()

# # ================================
# # Especially for Bales
# # ================================
# def get_row_1(pr_obj, x, h):
#     a, b, c, d = pr_obj.get_coef_tuple(x, h)
#     return a + b, c, d - a * h
#
#
# def get_row_n(pr_obj, x, h):
#     x -= h
#     a, b, c, d = pr_obj.get_coef_tuple(x, h)
#     return a, b - 2 * c / (h-2), d - c * 2 * h / (h - 2)
#
# p = lambda x: 0
# q = lambda x: -2 / (x**2 * (x + 1))
# f = lambda x: 0
#
# step = 0.1
# interval = [1, 2]
# pr = FiniteDifferenceMethod(interval, p, q, f)
# s = pr.solve(step, get_row_1(pr, interval[0], step), get_row_n(pr, interval[1], step),
#              left_value=lambda y, h: y + h,
#              right_value=lambda y, h: 2 * (h - y) / (h - 2))
# real_s = [(x[0], 1 / x[0] + 1) for x in s]
# print(real_s)
# print(s)
# # ================================
