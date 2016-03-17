import sys
import os.path
import math

my_libs_path = os.path.normpath(os.path.join(sys.path[0], "../"))
sys.path.append(my_libs_path)

try:
    from matrix import Matrix
except ImportError:
    my_libs_path = os.path.normpath(os.path.join(sys.path[0], "../../"))
    sys.path.append(my_libs_path)
    from matrix import Matrix


class JacobiEigenvalueAlgorithm:
    def __init__(self, mtrx):
        self.__mtrx = mtrx.copy()
        self.__need_logging = False

    def set(self, mtrx):
        self.__mtrx = mtrx.copy()

    def log(self, flag):
        self.__need_logging = flag

    @property
    def need_logging(self):
        return self.__need_logging

    @staticmethod
    def __find_max_above_diagonal(mtrx):
        max = math.fabs(mtrx[0][1])
        row, col = 0, 1
        for i in range(mtrx.size[0]):
            for j in range(i + 1, mtrx.size[1]):
                if math.fabs(mtrx[i][j]) > max:
                    max = math.fabs(mtrx[i][j])
                    row, col = i, j
        return row, col

    @staticmethod
    def __get_precision(mtrx):
        prec = 0
        for i in range(mtrx.size[0]):
            for j in range(i + 1, mtrx.size[1]):
                prec += mtrx[i][j] ** 2
        return math.sqrt(prec)

    @staticmethod
    def __get_rotation_matrix(mtrx, row, col):
        if mtrx[row][row] == mtrx[col][col]:
            angle = math.pi / 4
        else:
            angle = math.atan(2 * mtrx[row][col] / (mtrx[row][row] - mtrx[col][col]))
            angle /= 2
        s = math.sin(angle)
        c = math.cos(angle)
        u = Matrix.identity(mtrx.size)
        u[row][col] = -s
        u[col][row] = s
        u[row][row] = c
        u[col][col] = c
        return u

    def solve(self, precision):
        a_k = self.__mtrx.copy()
        u_k = None
        u = Matrix.identity(self.__mtrx.size)
        prec_func = JacobiEigenvalueAlgorithm.__get_precision
        iter_count = 0

        if self.need_logging:
            print("Solve with precision =", precision)
            print("A0 = ")
            print(a_k)
            print("precision =", prec_func(a_k), "\n")

        while prec_func(a_k) > precision:
            row, col = self.__find_max_above_diagonal(a_k)

            if self.need_logging:
                print("max element of A{0} is A[{1}][{2}] = {3}\n".format(iter_count, row+1, col+1, a_k[row][col]))

            u_k = JacobiEigenvalueAlgorithm.__get_rotation_matrix(a_k, row, col)
            u = u * u_k
            a_k = u_k.transpose() * a_k * u_k
            iter_count += 1

            if self.need_logging:
                print("U{} = ".format(iter_count - 1))
                print(u_k)
                print("A{} = ".format(iter_count))
                print(a_k)
                print("precision =", prec_func(a_k), "\n")

        if self.need_logging:
            print("U = ")
            print(u)
            print("Number of iterations =", iter_count, "\n")

        return a_k, u

    @staticmethod
    def split_to_vectors(mtrx):
        vectors = []
        for i in range(mtrx.size[1]):
            vectors.append(Matrix((mtrx.size[0], 1)))
            last = vectors[-1]
            for j in range(mtrx.size[0]):
                last[j][0] = mtrx[j][i]
        return vectors

    @staticmethod
    def split_to_values(mtrx):
        values = []
        for i in range(mtrx.size[0]):
            values.append(mtrx[i][i])
        return values

if __name__ == "__main__":
    pass
