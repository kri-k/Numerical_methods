import sys
import os.path
import math
import cmath

my_libs_path = os.path.normpath(os.path.join(sys.path[0], "../"))
sys.path.append(my_libs_path)

try:
    from matrix import Matrix
except ImportError:
    my_libs_path = os.path.normpath(os.path.join(sys.path[0], "../../"))
    sys.path.append(my_libs_path)
    from matrix import Matrix


class QREigenvalueAlgorithm:
    def __init__(self, mtrx):
        self.mtrx = mtrx.copy()
        self.__need_logging = False

    def set(self, mtrx):
        self.mtrx = mtrx.copy()

    @property
    def need_logging(self):
        return self.__need_logging

    def log(self, flag):
        self.__need_logging = flag

    @staticmethod
    def get_householder_transform(mtrx, col):
        sum_pow_2 = 0
        for i in range(col, mtrx.size[0]):
            sum_pow_2 += mtrx[i][col] ** 2
        v = Matrix((mtrx.size[0], 1))
        v[col][0] = mtrx[col][col] + math.copysign(math.sqrt(sum_pow_2), mtrx[col][col])
        for i in range(col + 1, mtrx.size[0]):
            v[i][0] = mtrx[i][col]
        h_mtrx = Matrix.identity(mtrx.size)
        transpose_v = v.transpose()
        tmp_mtrx = v * transpose_v
        tmp_mtrx *= 2 / (transpose_v * v)[0][0]
        h_mtrx -= tmp_mtrx
        return h_mtrx

    @staticmethod
    def QR_decomposition(mtrx):
        q = Matrix.identity(mtrx.size)
        a_k = mtrx.copy()
        for i in range(mtrx.size[1] - 1):
            h_k = QREigenvalueAlgorithm.get_householder_transform(a_k, i)
            q = q * h_k
            a_k = h_k * a_k
        return q, a_k  # R = a_k; QR = A (mtrx)

    @staticmethod
    def max_abs_in_col(mtrx, col, first_row=0):
        if not (0 <= col < mtrx.size[1] and 0 <= first_row < mtrx.size[0]):
            return 0
        max = math.fabs(mtrx[first_row][col])
        for i in range(first_row + 1, mtrx.size[0]):
            if math.fabs(mtrx[i][col]) > max:
                max = math.fabs(mtrx[i][col])
        return max

    @staticmethod
    def solve_quadratic_equation(a, b, c, d=0.):
        c -= d
        D = b**2 - 4*a*c
        tmp = 2*a
        x1 = (-b - cmath.sqrt(D)) / tmp
        x2 = (-b + cmath.sqrt(D)) / tmp
        return x1, x2

    def solve(self, precision):
        a_k = self.mtrx.copy()
        iter_count = 0
        n = self.mtrx.size[0]
        not_found = list(range(n - 1))
        need_delete = []
        eigenvalues = [None] * n
        while len(not_found) > 0:
            iter_count += 1
            q, r = QREigenvalueAlgorithm.QR_decomposition(a_k)
            a_k = r * q

            if self.need_logging:
                print("Q{} = ".format(iter_count))
                print(q)
                print("R{} = ".format(iter_count))
                print(r)
                print("A{0} = R{0} * Q{0} = ".format(iter_count))
                print(a_k)

            for col in not_found:
                if col in need_delete:
                    continue
                max = QREigenvalueAlgorithm.max_abs_in_col(a_k, col, first_row=col + 2)

                if self.need_logging:
                    print("Max |value| under(include) {} row in {} col = {}".format(col + 2, col, max))

                if max > precision:
                    continue
                if math.fabs(a_k[col + 1][col]) < precision:  # real number
                    eigenvalues[col] = a_k[col][col]
                    need_delete.append(col)

                    if self.need_logging:
                        print("Add l{} = {} as real eigenvalue".format(col, eigenvalues[col]))

                else:  # complex number
                    max = QREigenvalueAlgorithm.max_abs_in_col(a_k, col + 1, col + 2)
                    if max > precision or math.fabs(a_k[col + 1][col + 1]) < precision:
                        continue
                    j = col
                    x1, x2 = QREigenvalueAlgorithm.solve_quadratic_equation(1,
                                                                            -a_k[j][j] - a_k[j + 1][j + 1],
                                                                            a_k[j][j] * a_k[j + 1][j + 1],
                                                                            a_k[j][j + 1] * a_k[j + 1][j])
                    if eigenvalues[j] is not None and abs(x1 - eigenvalues[j]) < precision:
                        need_delete += [col, col + 1] if col < n-2 else [col]
                        if self.need_logging:
                            print("Add  l{} = {}, l{} = {} as complex eigenvalues\n".format(col, col+1, x1, x2))
                    elif self.need_logging:
                        if eigenvalues[j] is None:
                            print("First found complex eigenvalues l{} = {}, l{} = {}\n".format(col, x1, col+1, x2))
                        else:
                            print("Add next complex eigenvalues l{} = {}, l{} = {}\n".format(col, x1, col+1, x2))
                    eigenvalues[j] = x1
                    eigenvalues[j + 1] = x2
            for i in need_delete:
                not_found.remove(i)
            need_delete.clear()
        if eigenvalues[-1] is None or type(eigenvalues[-2]) is not complex:
            eigenvalues[-1] = a_k[n - 1][n - 1]
            if self.need_logging:
                print("Add l{} = {} as real eigenvalue".format(col, eigenvalues[col]))

        if self.need_logging:
            print("Number of iterations =", iter_count)

        return eigenvalues


if __name__ == "__main__":
    pass
    # a22 = 2.67
    # a33 = -2.59
    # a23 = 6.58
    # a32 = 2.79
    # print(QREigenvalueAlgorithm.solve_quadratic_equation(1,
    #                                                      -a22-a33,
    #                                                      a22*a33,
    #                                                      a23*a32))
