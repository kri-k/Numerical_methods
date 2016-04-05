from math import fabs
from libs.math_function import MathFunction
from libs.matrix import Matrix


class FixedPointIteration:
    def __init__(self, math_func):
        self.math_func = math_func

    def solve(self, precision, first_points, limiters, log=False):
        solutions = []
        for i in range(len(first_points)):
            solutions.append(self.__solve(precision, first_points[i], limiters[i], log))
        return solutions

    def __solve(self, prec, first_point, q, log):
        prev_x = first_point
        cur_x = self.math_func(prev_x)
        iter_count = 1
        while q / (1 - q) * fabs(cur_x - prev_x) > prec:
            prev_x = cur_x
            cur_x = self.math_func(cur_x)
            iter_count += 1
        if log:
            print("Number of iterations (with first point x0 = {}) = {}".format(first_point, iter_count))
        return cur_x


class SystemFixedPointIteration:
    def __init__(self, *functions):
        self.functions = functions

    def solve(self, precision, first_points, limiters, log=False):
        solutions = []
        for i in range(len(first_points)):
            solutions.append(self.__solve(precision, first_points[i], limiters[i], log))
        return solutions

    def __solve(self, prec, first_point, q, log):
        prev_x = Matrix((1, len(first_point)))
        prev_x.fill(*first_point)
        cur_x = self.__next_x(prev_x)
        iter_count = 1
        while q / (1 - q) * (cur_x - prev_x).transpose().norm_inf() > prec:
            prev_x = cur_x
            cur_x = self.__next_x(prev_x)
            iter_count += 1
        if log:
            print("Number of iterations (with first point x0 = {}) = {}".format(first_point, iter_count))
        return cur_x[0]

    def __next_x(self, cur_x):
        next_x = Matrix(cur_x.size)
        for i, f in enumerate(self.functions):
            next_x[0][i] = f(*cur_x[0])
        return next_x


if __name__ == "__main__":
    pass
