from inspect import getargspec
from libs.matrix import Matrix


# It will be cool to rewrite MathFunction with expressions tree
class MathFunction:
    def __init__(self, func=None, *derivates_first_order):
        if func is None:
            self.func = lambda x: x
            self.num_args = 1
        else:
            self.func = func
            self.num_args = len(getargspec(func).args)
        self.derivatives = \
            derivates_first_order[:min(self.num_args, len(derivates_first_order))]

    def __call__(self, *args):
        return self.func(*args)
        
    def derivative(self, var=0):
        if -1 < var < len(self.derivatives):
            return self.derivatives[var]
        else:
            return lambda *args: 0

    def get_derivatives(self):
        return self.derivatives

    @staticmethod
    def get_jacobian_matrix(functions, arg):
        n = len(functions)
        m = len(arg)
        jac = Matrix((n, m))
        for i in range(n):
            for j in range(m):
                jac[i][j] = functions[i].derivative(var=j)(*arg)
        return jac

if __name__ == "__main__":
    pass
