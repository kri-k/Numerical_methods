from libs.gauss_method import *


# Ordinary Least Squares
class OLS:
    def __init__(self, *points, order=1):
        self.points = list(points)
        self.coefficients = None
        self.order = order
        self.__update_coefficients()

    def add_points(self, *points):
        self.points += points
        self.__update_coefficients()

    def __update_coefficients(self):
        n = self.order
        N = len(self.points)
        mtrx = Matrix((n + 1, n + 1))
        vec = Matrix((n + 1, 1))
        x_list = [i[0] for i in self.points]
        x_sum_pows = [sum(map(lambda x: x**i, x_list)) for i in range(2*n + 1)]
        for k in range(n + 1):
            for i in range(n + 1):
                mtrx[k][i] = x_sum_pows[k + i]
            tmp = 0
            for i in range(N):
                tmp += self.points[i][1] * x_list[i]**k
            vec[k][0] = tmp
        self.coefficients = GaussMethod(mtrx, vec).solve().transpose()[0]

    def __call__(self, x):
        y = 0
        pow_x = 1
        for i in range(self.order + 1):
            y += self.coefficients[i] * pow_x
            pow_x *= x
        return y

    def error_sum_squares(self):
        sum = 0
        for p in self.points:
            sum += (self(p[0]) - p[1])**2
        return sum

if __name__ == "__main__":
    pass
