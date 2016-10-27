from lab_5_1.finite_difference.implicit_method import *
from lab_5_1.finite_difference.explicit_method import *


def crank_nicolson_fd(main_args,
                      boundary_left_args, boundary_right_args,
                      initial_func,
                      min_x, max_x,
                      max_t,
                      step_x, step_t,
                      boundary_approximation_func='first_order_two_points',
                      *args):
    """Solve parabolic partial differential equation with Crank - Nicolson finite difference method"""
    if len(args) == 2:
        explicit_solution = args[0]
        implicit_solution = args[1]
    else:
        l = [main_args,
             boundary_left_args, boundary_right_args,
             initial_func,
             min_x, max_x,
             max_t,
             step_x, step_t,
             boundary_approximation_func]
        explicit_solution = explicit_fd(*l)
        implicit_solution = implicit_fd(*l)
    crank_nicolson_solution = [list(map(lambda x, y: x / 2 + y / 2, explicit_solution[i], implicit_solution[i]))
                               for i in range(len(explicit_solution))]
    return crank_nicolson_solution


if __name__ == '__main__':
    pass
