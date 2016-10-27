from math import fabs, copysign
from functools import reduce
from libs.math_function import MathFunction
from libs.matrix import Matrix


class FixedPointIteration:
    def __init__(self, math_func):
        self.math_func = math_func

    def solve(self, precision, monotonicity_intervals, log=False):
        solutions = []
        for i in monotonicity_intervals:
            solutions.append(self.__solve(precision, i, log))
        return solutions

    def __solve(self, prec, monotonicity_interval, log):
        der = self.math_func.derivative()
        a = abs(der(monotonicity_interval[0]))
        b = abs(der(monotonicity_interval[1]))
        lam = copysign(1 / max(a, b), der(monotonicity_interval[0]))
        phi = lambda x: x - lam * self.math_func(x)
        q = abs(1 - min(a, b) / max(a, b))
        first_point = monotonicity_interval[0]
        prev_x = first_point
        cur_x = phi(prev_x)
        iter_count = 1
        while q / (1 - q) * fabs(cur_x - prev_x) > prec:
            prev_x = cur_x
            cur_x = phi(cur_x)
            iter_count += 1
        if log:
            print("Number of iterations (with first point x0 = {}) = {}".format(first_point, iter_count))
        return cur_x


__DEBUG = 0

def dpr(*args, **kwargs):
    if __DEBUG:
        print(*args, **kwargs)


class SystemFixedPointIteration:
    def __init__(self, *functions):
        self.functions = functions

    def solve(self, precision, left_boundaries, right_boundaries, log=False):
        solutions = []
        for i in range(len(left_boundaries)):
            solutions.append(self.__solve(precision,
                                          left_boundaries[i], right_boundaries[i],
                                          log))
        return solutions

    def __solve(self, prec, left_boundary, right_boundary, log):
        points = reduce((lambda ss, cs: (s + [c] for s in ss for c in cs)),
                        zip(left_boundary, right_boundary),
                        [[]])
        lam = -1
        min_jac = None
        first_point = None
        for p in points:
            tmp = MathFunction.get_jacobian_matrix(self.functions, p).norm_inf()
            if tmp > lam:
                lam = tmp
                # first_point = p
            if min_jac is None:
                min_jac = tmp
                first_point = p
            elif tmp < min_jac:
                min_jac = tmp
                first_point = p
        lam = 1 / lam
        dpr('lambda =', lam)
        dpr('min jacobian =', min_jac)
        dpr('=> first point =', first_point)

        prev_x = Matrix((1, len(first_point)))
        prev_x.fill(*first_point)
        # jac = MathFunction.get_jacobian_matrix(self.functions, monotonicity_interval[0])
        # lam = 1 / jac.norm_inf()
        phi_list = [lambda *x, f=self.functions[i], l=lam, i=i: x[i] - l*f(*x)
                    for i in range(len(self.functions))]
        q = 0.9
        # relax_coef = 0.1
        cur_x = self.__next_x(prev_x, phi_list)
        # prev_x *= 1 - relax_coef
        # cur_x *= relax_coef
        # cur_x += prev_x
        dpr('iter 1')
        dpr('prev_x =', prev_x[0])
        dpr('cur_x =', cur_x[0])
        iter_count = 1
        while q / (1 - q) * (cur_x - prev_x).transpose().norm_inf() > prec:
            dpr('prec = ', (cur_x - prev_x).transpose().norm_inf())
            prev_x = cur_x
            cur_x = self.__next_x(prev_x, phi_list)
            iter_count += 1
            dpr('next x = ', cur_x[0])
        if log:
            pass
            # print("Number of iterations (with first point x0 = {}) = {}".format(first_point, iter_count))
        return cur_x[0]

    def __next_x(self, cur_x, phi_list):
        next_x = Matrix(cur_x.size)
        for i, f in enumerate(phi_list):
            next_x[0][i] = f(*cur_x[0])
            # print('phi{}('.format(i + 1), cur_x[0], ') =', next_x[0][i])
        return next_x

if __name__ == "__main__":
    pass
