if __name__ == "__main__":
    from math import pi
    from libs.integration import *

    f = lambda x: 1 / (x**2 + 4)
    interval = (-2, 2)

    print('''f(x) = 1 / (x^2 + 4)
x in [-2, 2]
h1 = 1.0
h2 = 0.5''')

    print('''
================
Rectangle Method
----------------''')
    print('h = 1.0')
    print('F(x) =', RectangleMethod(f, interval).solve(1.0))
    print('h = 0.5')
    print('F(x) =', RectangleMethod(f, interval).solve(0.5))

    print('''
================
Trapezoidal Rule
----------------''')
    print('h = 1.0')
    print('F(x) =', TrapezoidalRule(f, interval).solve(1.0))
    print('h = 0.5')
    print('F(x) =', TrapezoidalRule(f, interval).solve(0.5))

    print('''
================
Simpson's Rule
----------------''')
    print('h = 1.0')
    print('F(x) =', SimpsonsRule(f, interval).solve(1.0))
    print('h = 0.5')
    print('F(x) =', SimpsonsRule(f, interval).solve(0.5))

    print('''
================
Runge Romberg Richardson Method
----------------''')
    print('Real value of F =', pi / 4, '\n')
    tmp = runge_romberg_richardson_method(RectangleMethod, f, interval, 1.0, 2)
    print('RRR method for Rectangle Method =', tmp)
    print('Rectangle Method Error =', abs(pi/4 - tmp), '\n')

    tmp = runge_romberg_richardson_method(TrapezoidalRule, f, interval, 1.0, 2)
    print('RRR method for Trapezoidal Rule Method =', tmp)
    print('Trapezoidal Rule Error =', abs(pi/4 - tmp), '\n')

    tmp = runge_romberg_richardson_method(SimpsonsRule, f, interval, 1.0, 2)
    print("RRR method for Simpson's Rule Method =", tmp)
    print("Simpson's Rule Error =", abs(pi/4 - tmp), '\n')

