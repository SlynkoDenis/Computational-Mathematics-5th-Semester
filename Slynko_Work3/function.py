import math


class Function:
    name = "sqrt(x) * exp(-x)"

    @staticmethod
    def get_max():
        return 0.5, math.sqrt(0.5) * math.exp(-0.5)

    @staticmethod
    def function(x):
        if x < 0:
            raise ValueError(f"{x} < 0")
        return math.sqrt(x) * math.exp(-x) - Function.get_max()[1] * 0.5

    @staticmethod
    def derivative(x):
        if x <= 0:
            raise ValueError(f"{x} < 0")
        return -1.0 * math.sqrt(x) * math.exp(-x) + 0.5 * math.exp(-x) / math.sqrt(x)

    @staticmethod
    def second_derivative(x):
        if x <= 0:
            raise ValueError(f"{x} < 0")
        return -1.0 * math.exp(-x) / math.sqrt(x) + math.sqrt(x) * math.exp(-x) - 0.25 * math.exp(-x) / math.sqrt(x) / x
