from libs.tridiagonal_matrix_algorithm import *


class SplineInterpolation:
    def __init__(self, *points):
        self.points = list(points)
        self.coefficients = None
        self.__update_coefficients()

    def add_points(self, *points):
        self.points += points
        self.__update_coefficients()

    def __update_coefficients(self):
        self.coefficients = []
        n = len(self.points) - 1
        h = [None] * (n + 1)
        for i in range(1, n + 1):
            h[i] = self.points[i][0] - self.points[i - 1][0]
        mtrx = Matrix((n - 1, 3))
        vec = Matrix((n - 1, 1))
        mtrx[0][1] = 2 * (h[1] + h[2])
        mtrx[0][2] = h[2]
        vec[0][0] = 3 * ((self.points[2][1] - self.points[1][1]) / h[2] -
                         (self.points[1][1] - self.points[0][1]) / h[1])
        for i in range(3, n):
            row = i - 2
            mtrx[row][0] = h[i - 1]
            mtrx[row][1] = 2 * (h[i - 1] + h[i])
            mtrx[row][2] = h[i]
            vec[row][0] = 3 * ((self.points[i][1] - self.points[i - 1][1]) / h[i] -
                               (self.points[i - 1][1] - self.points[i - 2][1]) / h[i - 1])
        row = n - 2
        mtrx[row][0] = h[n - 1]
        mtrx[row][1] = 2 * (h[n - 1] + h[n])
        vec[row][0] = 3 * ((self.points[n][1] - self.points[n - 1][1]) / h[n] -
                           (self.points[n - 1][1] - self.points[n - 2][1]) / h[n - 1])
        c_vec = TDMA(mtrx, vec).solve()
        c = [0.]
        for row in c_vec:
            c.append(row[0])
        for i in range(n - 1):
            a = self.points[i][1]
            b = (self.points[i + 1][1] - self.points[i][1]) / h[i + 1]
            b -= h[i + 1] * (c[i + 1] + 2*c[i]) / 3
            d = (c[i + 1] - c[i]) / (3 * h[i + 1])
            self.coefficients.append((a, b, c[i], d))
        a = self.points[n - 1][1]
        b = (self.points[n][1] - self.points[n - 1][1]) / h[n]
        b -= 2 * h[n] * c[n - 1] / 3
        d = -c[n - 1] / (3 * h[n])
        self.coefficients.append((a, b, c[-1], d))

    def __call__(self, x):
        interval_num = None
        for i in range(1, len(self.points)):
            if x <= self.points[i][0]:
                interval_num = i
                break
        if interval_num is None:
            raise ValueError('arg must be in interval ({}, {})'. format(self.points[0],
                                                                        self.points[-1]))
        interval_num -= 1
        a, b, c, d = self.coefficients[interval_num]
        sub_x = (x - self.points[interval_num][0])
        return a + \
            b * sub_x + \
            c * sub_x ** 2 + \
            d * sub_x ** 3

if __name__ == "__main__":
    pass
