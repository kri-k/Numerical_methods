from libs.derivative import Derivative

example = Derivative((0, 1),
                     (0.1, 1.1052),
                     (0.2, 1.2214),
                     (0.3, 1.3499),
                     (0.4, 1.4918))
print(example(0.2, order=1))
print(example(0.4, order=2))
