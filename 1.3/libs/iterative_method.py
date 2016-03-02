import sys
import os.path

my_libs_path = os.path.normpath(os.path.join(sys.path[0], "../"))
sys.path.append(my_libs_path)

from matrix import Matrix


class IterativeMethod:
    def __init__(self, mtrx, vec):
        self.__mtrx = mtrx
        self.__vector = vec
        self.alpha = None
        self.beta = None
        self.__need_logging = False

    def __str__(self):
        s = ""
        for i in range(self.__mtrx.size[0]):
            for j in range(self.__mtrx.size[1]):
                s += "{}*x{} + ".format(self.__mtrx[i][j], j)
            s = s[:-2] + "= {}\n".format(self.__vector[i][0])
        return s

    def get_alpha_mtrx(self):
        n = self.__mtrx.size[0]
        alpha = Matrix(self.__mtrx.size)
        for i in range(n):
            for j in range(n):
                if i != j:
                    alpha[i][j] = -self.__mtrx[i][j]/self.__mtrx[i][i]
                else:
                    alpha[i][j] = 0.
        return alpha

    def get_beta_vec(self):
        n = self.__vector.size[0]
        beta = Matrix(self.__vector.size)
        for i in range(n):
            beta[i][0] = self.__vector[i][0]/self.__mtrx[i][i]
        return beta

    @staticmethod
    def __get_precision_func(mtrx):
        m = mtrx.norm()
        if m > 1:
            def f(prev_vec, cur_vec):
                n = (cur_vec - prev_vec).norm()
                return m * n / (1 - m)
        else:
            def f(prev_vec, cur_vec):
                n = (cur_vec - prev_vec).norm()
                return n
        return f

    def solve(self, precision):
        if self.alpha is None:
            self.alpha = self.get_alpha_mtrx()
        if self.beta is None:
            self.beta = self.get_beta_vec()

        prec_func = self.__get_precision_func(self.alpha)

        if self.need_logging:
            print(self)
            print("Solve system with {} precision\n".format(precision))
            print("alpha matrix:", self.alpha, sep="\n")
            print("beta vector:", self.beta, sep="\n")

        prev_x = self.beta
        cur_x = self.beta + self.alpha * prev_x
        iter_count = 1

        if self.need_logging:
            print("X0:", prev_x, sep="\n")
            print("X1:", cur_x, sep="\n")
            print("precision = ", prec_func(prev_x, cur_x), end="\n\n")

        while precision < prec_func(prev_x, cur_x):
            iter_count += 1
            prev_x = cur_x
            cur_x = self.beta + self.alpha * cur_x
            if self.need_logging:
                print("X{}".format(iter_count), cur_x, sep="\n")
                print("precision = ", prec_func(prev_x, cur_x), end="\n\n")

        if self.need_logging:
            print("number of iterations:", iter_count, "\n")

        return cur_x

    def log(self, flag):
        self.__need_logging = flag

    @property
    def need_logging(self):
        return self.__need_logging


