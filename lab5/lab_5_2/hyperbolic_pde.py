"""Parabolic partial differential equation"""

# from lab_5_2.finite_difference import *


class HyperbolicPDE:
    """Solve hyperbolic partial differential equation

    Equation:
        d^2U/dt^2 = a^2 * d^2U/dx^2 + b * dU/dx + c * U + e * dU/dt + f(x, t)

        alpha * dU/dx(0, t) + beta * U(0, t) = phi_0(t)
        gamma * dU/dx(l, t) + delta * U(l, t) = phi_1(t)

        U(x, 0) = psi_1(x)
        dU/dt(x, 0) = psi_2(x)

        x in [0, l]
        t in [0, T]
    """

    def __init__(self,
                 a, b, c, f,
                 alpha, beta, phi_0,
                 gamma, delta, phi_1,
                 psi_1, psi_2,
                 min_x, max_x,
                 max_t):
        self.main_args = (a, b, c, e, f)
        self.boundary_left_args = (alpha, beta, phi_0)
        self.boundary_right_args = (gamma, delta, phi_1)
        self.initial_args = (psi_1, psi_2)
        self.max_x = max_x
        self.min_x = min_x
        self.max_t = max_t

    def solve(self, *args, **kwargs):
      pass
