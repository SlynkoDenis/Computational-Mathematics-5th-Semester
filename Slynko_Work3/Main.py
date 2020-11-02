import math
import function as function


class Main:
    f = function.Function()
    e = 0.001

    def __init__(self, epsilon=0.001):
        e = epsilon

    def newtons_method(self, start_approximation, func):
        while start_approximation > 0\
              and ((self.f.function(0)) * (self.f.function(start_approximation)) >= 0):
            start_approximation -= 0.005

        n = 0
        while abs(func.function(start_approximation)) > self.e:
            n += 1
            start_approximation -= (func.function(start_approximation)) / func.derivative(start_approximation)

        return start_approximation, n

    def fixed_point_iteration_method(self, start_approximation, func):
        while start_approximation > 0\
              and ((self.f.function(0)) * (self.f.function(start_approximation)) >= 0):
            start_approximation -= 0.005

        n = 0
        while abs(func.function(start_approximation)) > self.e:
            n += 1
            start_approximation += func.function(start_approximation)

        return start_approximation, n

    def fixed_point_iteration_method_reverse(self, start_approximation, func):
        while start_approximation > 0\
              and ((self.f.function(0)) * (self.f.function(start_approximation)) >= 0):
            start_approximation -= 0.005

        n = 0
        while abs(func.function(start_approximation)) > self.e:
            n += 1
            start_approximation -= func.function(start_approximation)

        return start_approximation, n

    def find_fwhm(self, method_left, method_right):
        max_x, max_value = function.Function.get_max()

        x_l, n_l = method_left(max_x - 0.4, self.f)
        print(f"Left roote is {x_l}")
        x_r, n_r = method_right(max_x + 0.2, self.f)
        print(f"Right roote is {x_r}")

        return x_r - x_l, n_l, n_r

    def show_results(self):
        print("Solving with Newton's method")
        r_n, n_l_n, n_r_n = self.find_fwhm(self.newtons_method, self.newtons_method)
        print("\nSolving with fixed-point iteration method")
        r_fpi, n_l_fpi, n_r_fpi = self.find_fwhm(self.fixed_point_iteration_method_reverse, self.fixed_point_iteration_method)
        print()

        print(f"For the fixed-point iteration method results are:\nx_r - x_l = {r_fpi}")
        print(f"Number of iterations for finding the left limit: {n_l_fpi}")
        print(f"Number of iterations for finding the right limit: {n_r_fpi}\n")

        print(f"For the Newton's method results are:\nx_r - x_l = {r_n}")
        print(f"Number of iterations for finding the left limit: {n_l_n}")
        print(f"Number of iterations for finding the right limit: {n_r_n}")


if __name__ == '__main__':
    main = Main()
    main.show_results()
