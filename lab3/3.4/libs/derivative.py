class Derivative:
    def __init__(self, *points):
        self.points = sorted(points, key=lambda x: x[0])
        self.first_ders = None
        self.second_ders = None
        self.__update_ders()

    def find_interval(self, x):
        for i, p in enumerate(self.points):
            if x <= p[0]:
                return i - 1

    def __update_ders(self):
        self.__update_first_ders()
        self.__update_second_ders()

    def __update_first_ders(self):
        self.first_ders = []
        for i in range(len(self.points) - 2):
            x0, y0 = self.points[i]
            x1, y1 = self.points[i + 1]
            x2, y2 = self.points[i + 2]
            c1 = (y1 - y0) / (x1 - x0)
            c2 = (y2 - y1) / (x2 - x1)
            c3 = x2 - x0
            c4 = -x0 - x1
            self.first_ders.append(lambda x, c1=c1, c2=c2, c3=c3, c4=c4: c1 + (c2 - c1) / c3 * (2 * x + c4))
        p0 = self.points[-2]
        p1 = self.points[-1]
        # at the right end we can calculate
        # first derivative only with first order of accuracy
        self.first_ders.append(lambda x, p0=p0, p1=p1: (p1[1] - p0[1]) / (p1[0] - p0[0]))

    def __update_second_ders(self):
        self.second_ders = []
        # at the right end we can't calculate
        # second derivative coz we have x_i & x_(i+1) but also needed x_(i+2)
        for i in range(len(self.points) - 2):
            x0, y0 = self.points[i]
            x1, y1 = self.points[i + 1]
            x2, y2 = self.points[i + 2]
            self.second_ders.append(((y2 - y1) / (x2 - x1) - (y1 - y0) / (x1 - x0)) / (x2 - x0))

    def first_order(self, x):
        i = self.find_interval(x)
        if i < len(self.points):
            return self.first_ders[i](x)
        raise ValueError('(((')

    def second_order(self, x):
        i = self.find_interval(x)
        if i < len(self.points) - 2:
            return self.second_ders[i]
        raise ValueError('(((')

    def __call__(self, x, order=1):
        if order == 1:
            return self.first_order(x)
        if order == 2:
            return self.second_order(x)
        return None
