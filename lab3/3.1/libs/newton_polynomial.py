class NewtonPolynomial:
    def __init__(self, *points):
        self.points = list(points)
        self.coefficients = [[] for _ in range(len(points))]
        self.__update_coefficients()

    def add_points(self, *points):
        self.points += points
        self.coefficients += [[] for _ in range(len(points))]
        self.__update_coefficients()

    def __update_coefficients(self):
        for i in range(len(self.coefficients)):
            self.__set_divided_differences(i)

    def __set_divided_differences(self, order):
        if order == 0:
            for i in range(len(self.coefficients[0]), len(self.points)):
                self.coefficients[0].append(self.points[i][1])
        else:
            for i in range(len(self.coefficients[order]), len(self.coefficients[order - 1]) - 1):
                tmp = self.coefficients[order - 1][i] - self.coefficients[order - 1][i + 1]
                tmp /= self.points[i][0] - self.points[i + order][0]
                self.coefficients[order].append(tmp)

    def __call__(self, x):
        mul = 1
        y = 0
        for i, p in enumerate(self.points):
            y += self.coefficients[i][0] * mul
            mul *= x - p[0]
        return y

if __name__ == "__main__":
    pass
