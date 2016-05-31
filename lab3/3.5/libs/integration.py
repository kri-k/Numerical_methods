class RectangleMethod:
    p = 2

    def __init__(self, f, interval):
        self.f = f
        self.interval = tuple(interval)

    def solve(self, step):
        solution = 0
        step_num = int((self.interval[1] - self.interval[0]) / step)
        prev_x = self.interval[0]
        for _ in range(step_num):
            cur_x = prev_x + step
            solution += self.f((cur_x + prev_x) / 2)
            prev_x = cur_x
        return solution * step


class TrapezoidalRule:
    p = 2

    def __init__(self, f, interval):
        self.f = f
        self.interval = tuple(interval)

    def solve(self, step):
        solution = 0
        step_num = int((self.interval[1] - self.interval[0]) / step)
        prev_x = self.interval[0]
        prev_y = self.f(prev_x)
        for _ in range(step_num):
            cur_x = prev_x + step
            cur_y = self.f(cur_x)
            solution += prev_y + cur_y
            prev_x = cur_x
            prev_y = cur_y
        return solution * step / 2


class SimpsonsRule:
    p = 4

    def __init__(self, f, interval):
        self.f = f
        self.interval = tuple(interval)

    def solve(self, step):
        solution = 0

        step_num = int((self.interval[1] - self.interval[0]) / (2 * step))
        prev_x = self.interval[0]
        prev_y = self.f(prev_x)
        for _ in range(step_num):
            solution += prev_y

            prev_x += step
            solution += 4 * self.f(prev_x)

            prev_x += step
            prev_y = self.f(prev_x)
            solution += prev_y
        return solution * step / 3


def runge_romberg_richardson_method(cls, f, interval, step, k):
    int_cls = cls(f, interval)
    sk = int_cls.solve(step)
    s = int_cls.solve(step / k)
    return s + (s - sk) / (k**cls.p - 1)


if __name__ == "__main__":
    pass
