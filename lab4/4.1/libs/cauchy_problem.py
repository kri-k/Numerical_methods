class CauchyProblem:
    def __init__(self, initial_value, *functions):
        self.func = tuple(functions)
        self.init_val = initial_value.copy()
        self.method_order = dict(
            euler=1,
            runge_kutta=4,
            adams=4)

    def solve(self, interval, step, method_name):
        switch = dict(
            euler=self.solve_euler,
            runge_kutta=self.solve_runge_kutta,
            adams=self.solve_adams)
        return switch[method_name](interval, step)

    def solve_euler(self, interval, step):
        sol = []
        xi = interval[0]
        yi = self.init_val.copy()
        sol.append((xi, yi))
        while xi + step <= interval[1]:
            yii = [yi[i] + step * self.func[i](xi, *yi) for i in range(len(yi))]
            xi += step
            sol.append((xi, yii))
            yi = yii[:]
        return sol

    def __K(self, xk, yk, step, i=0, prev_K=None):
        if i == 0:
            k_list = []
            l = [step * f(xk, *yk) for f in self.func]
            k_list.append(l)
            for j in range(1, 4):
                k_list.append(self.__K(xk, yk, step, j, prev_K=k_list[j-1]))
            return k_list
        elif i in (1, 2):
            new_yk = [yk[i] + 0.5 * prev_K[i] for i in range(len(yk))]
            l = [step * f(xk + 0.5 * step, *new_yk) for f in self.func]
            return l
        elif i == 3:
            new_yk = [yk[i] + prev_K[i] for i in range(len(yk))]
            l = [step * f(xk + step, *new_yk) for f in self.func]
            return l

    def solve_runge_kutta(self, interval, step):
        sol = []
        xi = interval[0]
        yi = self.init_val[:]
        sol.append((xi, yi))
        while xi + step <= interval[1]:
            k_list = self.__K(xi, yi, step)
            delta_yi = [(k_list[0][i] +
                         2*k_list[1][i] +
                         2*k_list[2][i] +
                         k_list[3][i]) / 6 for i in range(len(yi))]
            yii = [yi[i] + delta_yi[i] for i in range(len(yi))]
            xi += step
            sol.append((xi, yii))
            yi = yii[:]
        return sol

    def solve_adams(self, interval, step):
        sol = self.solve_runge_kutta((interval[0],
                                      interval[0] + step * 5), step)[:4]
        xi = sol[-1][0]
        yi = sol[-1][1]
        while xi + step <= interval[1]:
            yii = [yi[i] + step / 24 *
                   (55 * self.func[i](sol[-1][0], *sol[-1][1]) -
                    59 * self.func[i](sol[-2][0], *sol[-2][1]) +
                    37 * self.func[i](sol[-3][0], *sol[-3][1]) -
                    9 * self.func[i](sol[-4][0], *sol[-4][1])) for i in range(len(yi))]
            xi += step
            sol.append((xi, yii))
            yi = yii[:]
        return sol

    @staticmethod
    def runge_romberg(cauchy_problem_obj, interval, step, method_name):
        sol_2h = cauchy_problem_obj.solve(interval, step, method_name)
        sol_h = cauchy_problem_obj.solve((interval[0], interval[1] + step),
                                         step / 2,
                                         method_name)
        rr_sol = []
        d = 2 ** cauchy_problem_obj.method_order[method_name] - 1
        n = len(sol_h[0][1])
        for i in range(len(sol_2h)):
            rr_sol.append((sol_2h[i][0],
                          [sol_h[2*i][1][j] + (sol_h[2*i][1][j] - sol_2h[i][1][j]) / d for j in range(n)]))
        return rr_sol


if __name__ == "__main__":
    pass
