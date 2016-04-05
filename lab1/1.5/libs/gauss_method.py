import sys
import os.path

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

    def __str__(self):
        s = ""
        for i in range(self.__mtrx.size[0]):
            for j in range(self.__mtrx.size[1]):
                s += "{}*x{} + ".format(self.__mtrx[i][j], j+1)
            s = s[:-2] + "= {}\n".format(self.__vector[i][0])
        return s

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

        k = GaussMethod.__max_row_in_col(ak_mtrx, 0, 0)
        if k != 0:
            self.__permutations.append((0, k))
            ak_mtrx.exchange_rows(0, k)

        ak_mtrx = self.__update_L_mtrx(ak_mtrx, 0) * ak_mtrx

        for i in range(1, self.__mtrx.size[0]-1):
            k = GaussMethod.__max_row_in_col(ak_mtrx, i, i)
            if i != k:
                self.__permutations.append((i, k))
                ak_mtrx.exchange_rows(i, k)
                self.__L.exchange_rows(i, k)
                self.__L.exchange_cols(i, k)
            ak_mtrx = self.__update_L_mtrx(ak_mtrx, i) * ak_mtrx
        self.__U = ak_mtrx

    def solve(self):
        return self.__stage_2(self.__stage_1())

    def __stage_1(self):
        if self.__L is None or self.__U is None:
            self.__LU_decomposition()

        b = self.__vector.copy()

        for i in self.__permutations:
            b.exchange_rows(i[0], i[1])

        z = Matrix([self.__mtrx.size[0], 1])
        z[0][0] = b[0][0]

        for i in range(1, z.size[0]):
            sum = 0.
            for j in range(0, i):
                sum += self.__L[i][j] * z[j][0]
            z[i][0] = b[i][0] - sum
        return z

    def __stage_2(self, z):
        n = self.__mtrx.size[0]
        x = Matrix([n, 1])
        x[n-1][0] = z[n-1][0] / self.__U[n-1][n-1]

        for i in range(n-2, -1, -1):
            sum = 0
            for j in range(i+1, n):
                sum += self.__U[i][j] * x[j][0]
            x[i][0] = (z[i][0] - sum) / self.__U[i][i]
        return x

    @staticmethod
    def __max_row_in_col(mtrx, col, first_row):
        max_elem = abs(mtrx[first_row][col])
        row = first_row
        for i in range(first_row+1, mtrx.size[0]):
            if abs(mtrx[i][col]) > max_elem:
                max_elem = abs(mtrx[i][col])
                row = i
        return row

    @staticmethod
    def get_inverse_matrix(mtrx):
        gm = GaussMethod(mtrx)
        inv_mtrx = Matrix([mtrx.size[0], 0])
        b = Matrix([mtrx.size[0], 1])
        for i in range(mtrx.size[0]):
            b[i-1][0] = 0
            b[i][0] = 1
            gm.set_vector(b)
            inv_mtrx.append_col(gm.solve())
        return inv_mtrx

    @staticmethod
    def get_determinant(mtrx):
        gm = GaussMethod(mtrx)
        l, u = gm.get_LU()
        n = len(gm.__permutations)
        det = 1 if n % 2 == 0 else -1
        for i in range(u.size[0]):
            det *= u[i][i]
        return det
