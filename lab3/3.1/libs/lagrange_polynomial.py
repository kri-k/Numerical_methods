class LagrangePolynomial:
    def __init__(self, *points):
        self.points = list(points)
        self.coefficients = None
        self.__update_coefficients()

    def add_points(self, *points):
        self.points += points
        self.__update_coefficients()

    def __update_coefficients(self):
        self.coefficients = []
        for i in range(len(self.points)):
            tmp = 1
            for j in range(len(self.points)):
                if i == j:
                    continue
                tmp *= self.points[i][0] - self.points[j][0]
            self.coefficients.append(self.points[i][1] / tmp)

    def __mult(self, x):
        tmp = 1
        for xi in self.points:
            tmp *= x - xi[0]
        return tmp

    def __call__(self, x):
        y = 0
        for i, p in enumerate(self.points):
            y += self.coefficients[i] * self.__mult(x) / (x - p[0])
        return y

if __name__ == "__main__":
    pass
