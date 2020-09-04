import matplotlib.pyplot as plt
import math
import subprocess
import os
import sys
import functions

def numerical_diff_plus_h(function, h, x0):
    return (function(x0 + h) - function(x0)) / h

def numerical_diff_minus_h(function, h, x0):
    return (function(x0) - function(x0 - h)) / h

def numerical_diff_plus_minus_h(function, h, x0):
    return (function(x0 + h) - function(x0 - h)) / (2 * h)

def numerical_diff(function, h, x0):
    return 4.0 / 3.0 * numerical_diff_plus_minus_h(function, h, x0) -\
           1.0 / 3.0 * (function(x0 + 2 * h) -
                        function(x0 - 2 * h)) / (4 * h)

def numerical_diff_v2(function, h, x0):
    return 1.5 * numerical_diff_plus_minus_h(function, h, x0) -\
           0.6 * (function(x0 + 2 * h) -
                  function(x0 - 2 * h)) / (4 * h) + 0.1 * (function(x0 + 3 * h) -
                  function(x0 - 3 * h)) / (6 * h)

def make_and_save_plot(function_class, num_diff_functions, x0):
    if not os.path.isdir(os.path.join(os.path.dirname(sys.argv[0]), "Plots")):
        subprocess.check_output("mkdir Plots", shell=True)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for j in range(0, len(num_diff_functions)):
        x_list = [1 / math.pow(2, i) for i in range(0, 21)]
        y_list = [abs(function_class.derivative(x0) -
                      num_diff_functions[j](function_class.function, x_list[i], x0)) for i
                  in range(0, 21)]

        ax.plot(x_list, y_list, marker='o')

    ax.set_title(f"{function_class.name}", fontsize=20)
    ax.set_xlabel('step', fontsize=16)
    ax.set_ylabel('error', fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid()
    ax.legend(list(range(1, len(num_diff_functions) + 1)))
    fig.set_figheight(20)
    fig.set_figwidth(25)
    fig.savefig(os.path.join("Plots", f"{function_class.name}.png"))
