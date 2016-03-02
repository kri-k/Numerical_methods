import sys
import os.path
from .iterative_method import IterativeMethod

my_libs_path = os.path.normpath(os.path.join(sys.path[0], "../"))
sys.path.append(my_libs_path)

from matrix import Matrix


class SeidelMethod(IterativeMethod):
    def __init__(self, mtrx, vec):
        super(SeidelMethod, self).__init__(mtrx, vec)
        self.mtrx_C = None

    @staticmethod
    def __get_precision_func(mtrx):
        m = mtrx.norm()
        if m > 1:
            c = SeidelMethod.__get_upper_triangular_matrix(mtrx).norm()

            def f(prev_vec, cur_vec):
                n = (cur_vec - prev_vec).norm()
                return c * n / (1 - m)
        else:
            def f(prev_vec, cur_vec):
                n = (cur_vec - prev_vec).norm()
                return n
        return f

    @staticmethod
    def __get_upper_triangular_matrix(mtrx):
        m = Matrix(mtrx.size)
        for i in range(m.size[0]):
            for j in range(m.size[1]):
                if j >= i:
                    m[i][j] = mtrx[i][j]
        return m

    def solve(self, precision):
        if self.alpha is None:
            self.alpha = self.get_alpha_mtrx()
        if self.beta is None:
            self.beta = self.get_beta_vec()

        prec_func = SeidelMethod.__get_precision_func(self.alpha)

        if self.need_logging:
            print(self)
            print("Solve system with {} precision\n".format(precision))
            print("alpha matrix:", self.alpha, sep="\n")
            print("beta vector:", self.beta, sep="\n")

        prev_x = self.beta
        cur_x = self.__get_next_x(prev_x)
        iter_count = 1

        if self.need_logging:
            print("X0:", prev_x, sep="\n")
            print("X1:", cur_x, sep="\n")
            print("precision = ", prec_func(prev_x, cur_x), end="\n\n")

        while precision < prec_func(prev_x, cur_x):
            iter_count += 1
            prev_x = cur_x
            cur_x = self.__get_next_x(cur_x)
            if self.need_logging:
                print("X{}".format(iter_count), cur_x, sep="\n")
                print("precision = ", prec_func(prev_x, cur_x), end="\n\n")

        if self.need_logging:
            print("number of iterations:", iter_count, "\n")

        return cur_x

    def __get_next_x(self, cur_x):
        next_x = self.beta.copy()
        for i in range(next_x.size[0]):
            for j in range(i):
                next_x[i][0] += self.alpha[i][j] * next_x[j][0]
            for j in range(i, cur_x.size[0]):
                next_x[i][0] += self.alpha[i][j] * cur_x[j][0]
        return next_x
