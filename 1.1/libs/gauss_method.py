import sys
import os.path
from math import fabs

my_libs_path = os.path.normpath(os.path.join(sys.path[0], "../"))
sys.path.append(my_libs_path)

from matrix import Matrix


class GaussMethod:
    def __init__(self, mtrx=None, vec=None):
        self.__mtrx = mtrx.copy() if mtrx is not None else Matrix([1, 1])
        self.__vector = vec.copy() if vec is not None else Matrix([1, 1])
        self.__L = None
        self.__U = None
        self.__permutations = []
        self.__need_logging = False

    def __str__(self):
        s = ""
        for i in range(self.__mtrx.size[0]):
            for j in range(self.__mtrx.size[1]):
                s += "{}*x{} + ".format(self.__mtrx[i][j], j+1)
            s = s[:-2] + "= {}\n".format(self.__vector[i][0])
        return s

    def log(self, flag=False):
        self.__need_logging = flag

    def set_vector(self, vec):
        self.__vector = vec.copy()

    def set_matrix(self, mtrx):
        self.__mtrx = mtrx.copy()
        self.__LU_decomposition()

    def get_LU(self):
        if self.__L is None or self.__U is None:
            self.__LU_decomposition()
        return self.__L.copy(), self.__U.copy()

    def set_LU(self, L, U):
        self.__L = L.copy()
        self.__U = U.copy()

    @staticmethod
    def __get_mu(prev_a_mtrx, i, k):
        return prev_a_mtrx[i][k] / prev_a_mtrx[k][k]

    def __update_L_mtrx(self, prev_a_mtrx, k):
        mk = Matrix.identity(self.__mtrx.size)
        for i in range(k+1, mk.size[0]):
            mu = self.__get_mu(prev_a_mtrx, i, k)
            mk[i][k] = -mu
            self.__L[i][k] = mu
        return mk

    def __LU_decomposition(self):
        self.__L = Matrix.identity(self.__mtrx.size)
        ak_mtrx = self.__mtrx.copy()
        if self.need_logging:
            print("A0 = ")
            print(ak_mtrx)

        k = GaussMethod.__max_row_in_col(ak_mtrx, 0, 0)
        if k != 0:
            self.__permutations.append((0, k))
            ak_mtrx.exchange_rows(0, k)
            if self.need_logging:
                print("In A{} exchange {} and {} rows\n".format(0, 0, k))

        if self.need_logging:
            print("PA{} = ".format(0))
            print(ak_mtrx)
            print("L{} = ".format(0))
            print(self.__L)

        ak_mtrx = self.__update_L_mtrx(ak_mtrx, 0) * ak_mtrx

        if self.need_logging:
            print("A{} = ".format(1))
            print(ak_mtrx)
            print("L{} = ".format(1))
            print(self.__L)

        for i in range(1, self.__mtrx.size[0]-1):
            k = GaussMethod.__max_row_in_col(ak_mtrx, i, i)
            if i != k:
                self.__permutations.append((i, k))
                ak_mtrx.exchange_rows(i, k)
                self.__L.exchange_rows(i, k)
                self.__L.exchange_cols(i, k)
                if self.need_logging:
                    print("In A{} exchange {} and {} rows".format(i, i, k))
                    print("In L{} exchange {} and {} rows".format(i, i, k))
                    print("In L{} exchange {} and {} columns\n".format(i, i, k))

            if self.need_logging:
                print("PA{} = ".format(i))
                print(ak_mtrx)
                print("L{} = ".format(i))
                print(self.__L)

            ak_mtrx = self.__update_L_mtrx(ak_mtrx, i) * ak_mtrx

            if self.need_logging:
                print("A{} = ".format(i+1))
                print(ak_mtrx)
                print("L{} = ".format(i+1))
                print(self.__L)

        if self.need_logging:
            print("U = A{}\n".format(self.__mtrx.size[0]-1))

        self.__U = ak_mtrx

        if self.need_logging:
            print("L * U = ")
            print(self.__L * self.__U)
            pa = self.__mtrx.copy()
            for i in self.__permutations:
                pa.exchange_rows(i[0], i[1])
            print("PA = ")
            print(pa)

    def solve(self):
        if self.need_logging:
            print(self)
            print("Solve Ax = b")
            print("Transform to PAx = Pb => LUx = Pb\n")
        return self.__stage_2(self.__stage_1())

    def __stage_1(self):
        if self.__L is None or self.__U is None:
            self.__LU_decomposition()

        if self.need_logging:
            print("Solve Lz = Pb\n")

        b = self.__vector.copy()

        if self.need_logging:
            print("b = ")
            print(b)

        for i in self.__permutations:
            b.exchange_rows(i[0], i[1])

        if self.need_logging:
            print("Pb = ")
            print(b)

        z = Matrix([self.__mtrx.size[0], 1])
        z[0][0] = b[0][0]

        if self.need_logging:
            print("z[0] = b[0] =", b[0][0], "\n")

        for i in range(1, z.size[0]):
            sum = 0.
            for j in range(0, i):
                sum += self.__L[i][j] * z[j][0]
            z[i][0] = b[i][0] - sum

            if self.need_logging:
                str_theoretically = ''
                str_practically = ''
                for j in range(0, i):
                    str_theoretically += "L[{}][{}]*z[{}] + ".format(i, j, j)
                    str_practically += "{}*{} + ".format(self.__L[i][j], z[j][0])
                str_theoretically = "z[{}] = b[{}] - ".format(i, i) + str_theoretically[:-3]
                str_practically = "z[{}] = {} - ".format(i, b[i][0]) + str_practically[:-3]
                print(str_theoretically)
                print(str_practically)
                print("z[{}] = {}\n".format(i, z[i][0]))

        return z

    def __stage_2(self, z):
        if self.need_logging:
            print("Solve Ux = z\n")

        n = self.__mtrx.size[0]
        x = Matrix([n, 1])
        x[n-1][0] = z[n-1][0] / self.__U[n-1][n-1]

        if self.need_logging:
            print("x[{0}] = z[{0}] / U[{0}][{0}] = {1} / {2} = {3}\n".
                  format(n-1, z[n-1][0], self.__U[n-1][n-1], x[n-1][0]))

        for i in range(n-2, -1, -1):
            sum = 0
            for j in range(i+1, n):
                sum += self.__U[i][j] * x[j][0]
            x[i][0] = (z[i][0] - sum) / self.__U[i][i]

            if self.need_logging:
                str_theoretically = ''
                str_practically = ''
                for j in range(i+1, n):
                    str_theoretically += "U[{}][{}]*x[{}] + ".format(i, j, j)
                    str_practically += "{}*{} + ".format(self.__U[i][j], x[j][0])
                str_theoretically = "x[{0}] = (z[{0}] - ".format(i) + str_theoretically[:-3] + \
                    ") / U[{0}][{0}]".format(i)
                str_practically = "x[{}] = ({} - {}) / {}".format(i, z[i][0], sum, self.__U[i][i])
                print(str_theoretically)
                print(str_practically)
                print("x[{}] = {}\n".format(i, x[i][0]))

        return x

    @staticmethod
    def __max_row_in_col(mtrx, col, first_row):
        max_elem = fabs(mtrx[first_row][col])
        row = first_row
        for i in range(first_row+1, mtrx.size[0]):
            if fabs(mtrx[i][col]) > max_elem:
                max_elem = fabs(mtrx[i][col])
                row = i
        return row

    @staticmethod
    def get_inverse_matrix(mtrx, log=False):
        gm = GaussMethod(mtrx)
        gm.log(log)
        inv_mtrx = Matrix([mtrx.size[0], 0])
        b = Matrix([mtrx.size[0], 1])
        for i in range(mtrx.size[0]):
            b[i-1][0] = 0
            b[i][0] = 1
            gm.set_vector(b)
            inv_mtrx.append_col(gm.solve())
        if log:
            print("Check the correctness of the answer: ")
            print("A * A^(-1) = ")
            print(mtrx * inv_mtrx)
            print()
        return inv_mtrx

    @staticmethod
    def get_determinant(mtrx, log=False):
        gm = GaussMethod(mtrx)
        gm.log(log)
        l, u = gm.get_LU()
        n = len(gm.__permutations)
        det = 1 if n % 2 == 0 else -1
        for i in range(u.size[0]):
            det *= u[i][i]
        return det

    @property
    def need_logging(self):
        return self.__need_logging

