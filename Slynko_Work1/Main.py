import math
import functions
import plot

class Main:
    functions_list = None
    derivative_methods = None

    def __init__(self):
        self.functions_list = [functions.SinX2, functions.CosSinX, functions.ExpSinCosX,
                               functions.LnXPlus3, functions.SqrtXPlus3]
        self.derivative_methods = [plot.numerical_diff_plus_h, plot.numerical_diff_minus_h,
                                   plot.numerical_diff_plus_minus_h, plot.numerical_diff,
                                   plot.numerical_diff_v2]

    def make_all_plots(self):
        for i in range(0, len(self.functions_list)):
            plot.make_and_save_plot(self.functions_list[i], self.derivative_methods, 5.0)


if __name__ == '__main__':
    main = Main()
    main.make_all_plots()