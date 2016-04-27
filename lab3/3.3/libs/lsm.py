from libs.gauss_method import *


# Least Squares Method
class LSM:
    def __init__(self, *points, order=1, basis_funcs=None):
        self.points = list(points)
        self.coefficients = None
        self.order = order
        if basis_funcs is None:
            self.basis = [lambda x, i=i: x**i for i in range(order + 1)]
        else:
            self.basis = basis_funcs.copy()
        self.__update_coefficients()

    def add_points(self, *points):
        self.points += points
        self.__update_coefficients()

    def __update_coefficients(self):
        n = len(self.basis)
        N = len(self.points)
        mtrx = Matrix((N, n))
        vec = Matrix((N, 1))
        vec.fill(*[p[1] for p in self.points])
        mtrx_elems = []
        for p in self.points:
            mtrx_elems += [f(p[0]) for f in self.basis]
        mtrx.fill(*mtrx_elems)
        mtrx_t = mtrx.transpose()
        mtrx = mtrx_t * mtrx
        vec = mtrx_t * vec
        self.coefficients = GaussMethod(mtrx, vec).solve().transpose()[0]

    def __call__(self, x):
        y = 0
        for i in range(len(self.basis)):
            y += self.coefficients[i] * self.basis[i](x)
        return y

    def error_sum_squares(self):
        sum = 0
        for p in self.points:
            sum += (self(p[0]) - p[1])**2
        return sum

if __name__ == "__main__":
    pass
