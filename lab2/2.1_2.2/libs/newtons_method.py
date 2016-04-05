from math import fabs
from libs.math_function import MathFunction
from libs.matrix import Matrix
from libs.gauss_method import GaussMethod


class NewtonsMethod:
    def __init__(self, math_func):
        self.math_func = math_func

    def solve(self, precision, first_points, log=False):
        solutions = []
        for fp in first_points:
            solutions.append(self.__solve(precision, fp, log))
        return solutions

    def __solve(self, prec, first_point, log):
        prev_x = first_point
        cur_x = prev_x - self.math_func(prev_x) / self.math_func.derivative()(prev_x)
        iter_count = 1
        while fabs(cur_x - prev_x) > prec:
            prev_x = cur_x
            cur_x = prev_x - self.math_func(prev_x) / self.math_func.derivative()(prev_x)
            iter_count += 1
        if log:
            print("Number of iterations (with first point x0 = {}) = {}".format(first_point, iter_count))
        return cur_x


class SystemNewtonsMethod:
    def __init__(self, *functions):
        self.functions = functions

    def solve(self, precision, first_points, log=False):
        solutions = []
        for fp in first_points:
            solutions.append(self.__solve(precision, fp, log))
        return solutions

    def __solve(self, prec, first_point, log):
        prev_x = Matrix((1, len(first_point)))
        prev_x.fill(*first_point)
        cur_x = self.get_delta(prev_x).transpose() + prev_x
        iter_count = 1
        while (cur_x - prev_x).norm_inf() > prec:
            prev_x = cur_x
            cur_x = self.get_delta(prev_x).transpose() + prev_x
            iter_count += 1
        if log:
            print("Number of iterations (with first point x0 = {}) = {}".format(first_point, iter_count))
        return cur_x[0]

    def get_delta(self, cur_x):
        n = cur_x.size[1]
        a = Matrix((n, n))
        b = Matrix((n, 1))
        elems_a = []
        elems_b = []
        for f in self.functions:
            elems_b.append(-f(*cur_x[0]))
            der = f.get_derivatives()
            for df in der:
                elems_a.append(df(*cur_x[0]))
        a.fill(*elems_a)
        b.fill(*elems_b)
        return GaussMethod(a, b).solve()

if __name__ == "__main__":
    pass
