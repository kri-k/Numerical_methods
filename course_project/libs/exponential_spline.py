import math
from libs.tridiagonal_matrix_algorithm import *


__DEBUG = 1

def dpr(*args, **kwargs):
    if __DEBUG:
        print(*args, **kwargs)


class ExpSpline:
    def __init__(self, *points):
        self.points = list(points)
        self.tension_parameters = [0.0001 for _ in range(len(points) - 1)]

        self.h = self.__get_interval_dists(points)
        self.left_border_derivative = 0
        self.right_border_derivative = 0

        self.__set_tension_parameters()
        self.coefficients = self.__get_coefficients()

    def add_points(self, *points):
        self.points += points
        self.h += self.__get_interval_dists(points)
        self.tension_parameters = [0 for _ in range(len(points) - 1)]

        self.__set_tension_parameters()
        self.coefficients = self.__get_coefficients()

    @staticmethod
    def __get_interval_dists(points):
        return [points[i][0] - points[i-1][0] for i in range(1, len(points))]

    def __get_mtrx(self):
        n = len(self.points)
        mtrx = Matrix((n, 3))
        mtrx[0][1] = self.__d(0)
        mtrx[0][2] = self.__e(0)
        for i in range(1, n - 1):
            mtrx[i][0] = self.__e(i - 1)
            mtrx[i][1] = self.__d(i - 1) + self.__d(i)
            mtrx[i][2] = self.__e(i)
        mtrx[n-1][0] = self.__e(n-2)
        mtrx[n-1][1] = self.__d(n-2)
        return mtrx

    def __get_vec(self):
        n = len(self.points)
        vec = Matrix((n, 1))
        vec.fill(*[self.__b(i) for i in range(n)])
        return vec

    def __d(self, i):
        p = self.tension_parameters[i]
        h = self.h[i]
        s = math.sinh(p * h)
        c = math.cosh(p * h)
        return (p * c / s - 1 / h) / p**2

    def __e(self, i):
        p = self.tension_parameters[i]
        h = self.h[i]
        s = math.sinh(p * h)
        return (1 / h - p / s) / p**2

    def __b(self, i):
        if i == 0:
            return (self.points[1][1] - self.points[0][1]) / self.h[0] - \
                   self.left_border_derivative
        if i == len(self.points) - 1:
            return self.right_border_derivative - \
                   (self.points[-1][1] - self.points[-2][1]) / self.h[-1]
        b = (self.points[i + 1][1] - self.points[i][1]) / self.h[i]
        b -= (self.points[i][1] - self.points[i - 1][1]) / self.h[i - 1]
        return b

    def __get_coefficients(self):
        m = self.__get_mtrx()
        v = self.__get_vec()
        return TDMA(m, v).solve().transpose()[0]

    def __lam(self, i, coef):
        if i == 0:
            return max(abs(self.__b(0)),
                       self.__d(0) * abs(coef[0])) / coef[1]
        if i == len(self.points) - 1:
            n = len(self.points) - 1
            return max(abs(self.__b(n)),
                       self.__d(n - 1) * abs(coef[n])) / coef[n - 1]
        tmp = max(abs(self.__b(i)), (self.__d(i - 1) + self.__d(i)) * abs(coef[i]))
        tmp /= 2 * max(abs(coef[i - 1]), abs(coef[i + 1]))
        return tmp

    def __ph(self, i, lam):
        return max(abs(lam * self.h[i]) ** -0.5, self.tension_parameters[i])

    def __set_tension_parameters(self):
        relax_coef = 0.01
        max_delta = 20
        max_iter_num = max_delta // relax_coef

        coef = self.__get_coefficients()

        b = [self.__b(i) for i in range(len(coef))]

        dpr('===================')
        dpr('b =', ' '.join(map(lambda x: str(x)[:7] if x >= 0 else str(x)[:8], b)))
        dpr("t'' =", ' '.join(map(lambda x: str(x)[:7] if x >= 0 else str(x)[:8], coef)))

        for i in range(0, len(coef)):
            dpr('===================')
            dpr(i)
            dpr('b =', ' '.join(map(lambda x: str(x)[:7] if x >= 0 else str(x)[:8], b)))
            dpr("t'' =", ' '.join(map(lambda x: str(x)[:7] if x >= 0 else str(x)[:8], coef)))

            s = 0
            if b[i] * b[i-1] <= 0:
                dpr('\tb[{}] * b[{}] = {} <= 0: CONTINUE'.format(i, i-1, b[i] * b[i-1]))
                continue

            if coef[i] * b[i] < 0:
                dpr("\tt''[{0}] * b[{0}] = {1} < 0: ITERATE".format(i, coef[i] * b[i]))

            while coef[i] * b[i] < 0:
                # break
                s += 1

                lam = self.__lam(i, coef)
                self.tension_parameters[i - 1] += \
                    relax_coef * (self.__ph(i - 1, lam) - self.tension_parameters[i - 1])
                self.tension_parameters[i] += \
                    relax_coef * (self.__ph(i, lam) - self.tension_parameters[i])

                coef = self.__get_coefficients()
                if s > max_iter_num:
                    dpr('OVERDOZE!!!')
                    break

            dpr("\t\tt''[{0}] * b[{0}] = {1} >= 0: {2}".format(i, coef[i] * b[i], coef[i] * b[i] >= 0))
        dpr('p = ', ' '.join(map(lambda x: str(x)[:7] if x >= 0 else str(x)[:8], self.tension_parameters)))

    def __call__(self, x):
        n = None
        for i in range(1, len(self.points)):
            if x <= self.points[i][0]:
                n = i
                break
        if n is None:
            raise ValueError('arg must be in interval ({}, {})'. format(self.points[0],
                                                                        self.points[-1]))
        n -= 1
        sub_x_up = self.points[n+1][0] - x
        sub_x_dn = x - self.points[n][0]
        c = self.coefficients
        p = self.tension_parameters[n]
        ans = c[n] * math.sinh(p * sub_x_up)
        ans += c[n+1] * math.sinh(p * sub_x_dn)
        ans /= p ** 2 * math.sinh(p * self.h[n])
        ans += (self.points[n][1] - c[n] / p**2) * sub_x_up / self.h[n]
        ans += (self.points[n+1][1] - c[n+1] / p**2) * sub_x_dn / self.h[n]
        return ans

if __name__ == "__main__":
    pass
