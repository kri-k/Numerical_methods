from libs.tridiagonal_matrix_algorithm import *
from libs.cauchy_problem import CauchyProblem

__DEBUG = 0


def dpr(*args, **kwargs):
    if __DEBUG:
        print(*args, **kwargs)


class ShootingMethod:
    def __init__(self, interval, initial_value_funcs, functions, right_boundary_func):
        self.interval = tuple(interval)
        self.init_val_f = tuple(initial_value_funcs)
        self.func = tuple(functions)
        self.right_f = right_boundary_func

    def solve(self, step, first_params_pair, precision):
        if len(first_params_pair) != 2:
            return None

        interval = list(self.interval)
        # to be sure we found value at the right side
        # (this may be wrong due to errors of calculation)
        interval[1] += step / 2

        init_vals = [[f(param) for f in self.init_val_f] for param in first_params_pair]
        prev_param = first_params_pair[0]
        cur_param = first_params_pair[1]

        prev_right_val = self.right_f(CauchyProblem(init_vals[0], *self.func).solve_runge_kutta(interval, step)[-1][1])
        cur_sol = CauchyProblem(init_vals[1], *self.func).solve_runge_kutta(interval, step)
        cur_right_val = self.right_f(cur_sol[-1][1])

        dpr('F = %s' % prev_right_val)
        dpr('F = %s' % cur_right_val)

        while abs(cur_right_val) > precision:
            try:
                tmp = cur_param - (cur_param - prev_param) * cur_right_val / (cur_right_val - prev_right_val)
            except ZeroDivisionError as e:
                print(e)
                print(cur_right_val, prev_right_val)
                return cur_sol

            dpr('{0} - ({0} - {1}) * {2} / ({2} - {3})'.format(cur_param, prev_param, cur_right_val, prev_right_val))
            dpr('n = %s' % tmp)

            new_init_vals = [f(tmp) for f in self.init_val_f]
            dpr('init_vals = ', new_init_vals)

            prev_param = cur_param
            prev_right_val = cur_right_val

            cur_param = tmp
            cur_sol = CauchyProblem(new_init_vals, *self.func).solve_runge_kutta(interval, step)
            dpr('cur right sol = ', cur_sol[-1])
            cur_right_val = self.right_f(cur_sol[-1][1])

            dpr('F = %s' % cur_right_val)

        return cur_sol

    def runge_romberg(self, step, first_params_pair, precision):
        sol_2h = self.solve(step, first_params_pair, precision)
        sol_h = self.solve(step / 2, first_params_pair, precision)
        rr_sol = []
        n = len(sol_h[0][1])
        d = 2 ** 4 - 1
        for i in range(len(sol_2h)):
            rr_sol.append((sol_2h[i][0],
                          [sol_h[2*i][1][j] + (sol_h[2*i][1][j] - sol_2h[i][1][j]) / d for j in range(n)]))
        return rr_sol


class FiniteDifferenceMethod:
    def __init__(self, interval, p_x, q_x, f_x):
        self.interval = interval
        self.f = (p_x, q_x, f_x)

    def solve(self, step, row_1, row_n, left_value, right_value):
        h = step
        n = int((self.interval[1] - self.interval[0]) / step) + 1
        m = Matrix((n - 2, 3))
        v = Matrix((n - 2, 1))
        x_k = [self.interval[0] + step * i for i in range(n)]
        m[0][1] = row_1[0]
        m[0][2] = row_1[1]
        v[0][0] = row_1[2]
        for k in range(1, n-3):
            x = x_k[k]
            m[k][0] = self.A(x, h)
            m[k][1] = self.B(x, h)
            m[k][2] = self.C(x, h)
            v[k][0] = self.D(x, h)
        m[-1][0] = row_n[0]
        m[-1][1] = row_n[1]
        v[-1][0] = row_n[2]
        sol = TDMA(m, v).solve().transpose()[0]
        sol.insert(0, left_value(sol[0], h))
        sol.append(right_value(sol[-1], h))
        return list(zip(x_k, sol))

    def A(self, x, h):
        return 1 - self.f[0](x) * h / 2

    def B(self, x, h):
        return -2 + h**2 * self.f[1](x)

    def C(self, x, h):
        return 1 + self.f[0](x) * h / 2

    def D(self, x, h):
        return h**2 * self.f[2](x)

    def get_coef_tuple(self, x, h):
        return self.A(x, h), self.B(x, h), self.C(x, h), self.D(x, h)

    def runge_romberg(self, step, row_1, row_n, left_value, right_value):
        sol_2h = self.solve(step,
                            row_1(self, self.interval[0], step),
                            row_n(self, self.interval[1], step),
                            left_value,
                            right_value)
        sol_h = self.solve(step / 2,
                           row_1(self, self.interval[0], step / 2),
                           row_n(self, self.interval[1], step / 2),
                           left_value,
                           right_value)
        rr_sol = []
        d = 2 ** 1 - 1
        for i in range(len(sol_2h)):
            rr_sol.append((sol_2h[i][0],
                          sol_h[2*i][1] + (sol_h[2*i][1] - sol_2h[i][1]) / d))
        return rr_sol


if __name__ == "__main__":
    pass
