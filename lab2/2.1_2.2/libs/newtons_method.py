from math import fabs
from libs.math_function import MathFunction
from libs.matrix import Matrix
from libs.gauss_method import GaussMethod


class NewtonsMethod:
    def __init__(self, math_func):
        self.math_func = math_func

    def solve(self, precision, monotonicity_intervals, log=False):
        solutions = []
        for i in monotonicity_intervals:
            solutions.append(self.__solve(precision, i, log))
        return solutions

    def __solve(self, prec, monotonicity_interval, log):
        first_point = (monotonicity_interval[0] + monotonicity_interval[1]) / 2
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

    def solve(self, precision, monotonicity_intervals, log=False):
        solutions = []
        for i in monotonicity_intervals:
            solutions.append(self.__solve(precision, i, log))
        return solutions

    def __solve(self, prec, monotonicity_interval, log):
        """
         monotonicity_interval: two points in Rn
        """
        first_point = list(map(lambda x, y: (x + y) / 2, *monotonicity_interval))
        prev_x = Matrix((1, len(first_point)))
        prev_x.fill(*first_point)
        cur_x = self.get_delta(prev_x) + prev_x
        iter_count = 1
        while (cur_x - prev_x).norm_inf() > prec:
            prev_x = cur_x
            cur_x = self.get_delta(prev_x) + prev_x
            iter_count += 1
        if log:
            print("Number of iterations (with first point x0 = {}) = {}".format(first_point, iter_count))
        return cur_x[0]

    def get_delta(self, cur_x):
        n = cur_x.size[1]
        a = MathFunction.get_jacobian_matrix(self.functions, cur_x[0])
        b = Matrix((n, 1))
        elems_b = []
        for f in self.functions:
            elems_b.append(-f(*cur_x[0]))
        b.fill(*elems_b)
        return GaussMethod(a, b).solve().transpose()

if __name__ == "__main__":
    pass
