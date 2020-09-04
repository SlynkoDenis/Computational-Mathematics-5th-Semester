import math

class SinX2:
    name = "sin(x^2)"

    @staticmethod
    def function(x):
        return math.sin(x * x)

    @staticmethod
    def derivative(x):
        return math.cos(x * x) * 2.0 * x

class CosSinX:
    name = "cos(sin(x))"

    @staticmethod
    def function(x):
        return math.cos(math.sin(x))

    @staticmethod
    def derivative(x):
        return -1.0 * math.sin(math.sin(x)) * math.cos(x)

class ExpSinCosX:
    name = "exp(sin(cos(x)))"

    @staticmethod
    def function(x):
        return math.exp(math.sin(math.cos(x)))

    @staticmethod
    def derivative(x):
        return -1.0 * ExpSinCosX.function(x) * math.cos(math.cos(x)) * math.sin(x)

class LnXPlus3:
    name = "ln(x+3)"

    @staticmethod
    def function(x):
        return math.log(x + 3.0)

    @staticmethod
    def derivative(x):
        return 1.0 / (x + 3.0)

class SqrtXPlus3:
    name = "sqrt(x+3)"

    @staticmethod
    def function(x):
        return math.sqrt(x + 3.0)

    @staticmethod
    def derivative(x):
        return 0.5 / SqrtXPlus3.function(x)
