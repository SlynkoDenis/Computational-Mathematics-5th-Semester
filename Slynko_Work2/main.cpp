#include <iostream>
#include "matrix.h"
#include <chrono>

int main() {
    int dim = 0;
    std::cin >> dim;
    if (dim <= 0) {
        throw std::invalid_argument("dimension must be positive");
    }
    if (dim < 9) {
        throw std::invalid_argument("dimension doesn't suit the particular problem");
    }

    Matrix<long double> A(dim, dim);
    for (size_t i = 1; i <= 4; ++i) {
        for (size_t j = 1; j <= i + 4; ++j) {
            if (i == j) {
                A.at(i, i) = 20.0L;
                A.at(dim + 1 - i, dim + 1 - i) = 20.0L;
            } else {
                A.at(i, j) = 1.0L;
                A.at(dim + 1 - i, dim + 1 - j) = 1.0L;
            }
        }
    }
    for (size_t i = 5; i <= dim - 4; ++i) {
        for (size_t j = i - 4; j <= i + 4; ++j) {
            if (i == j) {
                A.at(i, i) = 20.0L;
            } else {
                A.at(i, j) = 1.0L;
            }
        }
    }

    Matrix<long double> f(dim, 1);
    for (size_t i = 1; i <= dim; ++i) {
        f.at(i, 1) = static_cast<long double>(i);
    }

    // Condition number of the matrix A
    long double mu = condition_number(A);

    std::chrono::steady_clock::time_point begin_time = std::chrono::steady_clock::now();

    // Gaussian elimination solution
    decltype(auto) gaussian_solution = gaussian_elimination(A, f);
    decltype(auto) gaussian_residual = f - A * gaussian_solution;

    std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
    std::cout << "Time elapsed for Gaussian elimination = " <<
                 std::chrono::duration_cast<std::chrono::microseconds>(end_time - begin_time).count() << "[µs]" << std::endl;

    begin_time = std::chrono::steady_clock::now();

    // Gauss-Seidel method solution
    decltype(auto) seidel_solution = gauss_seidel_method(A, f);
    decltype(auto) seidel_residual = f - A * seidel_solution;

    end_time = std::chrono::steady_clock::now();
    std::cout << "Time elapsed for Gauss-Seidel method = " <<
                 std::chrono::duration_cast<std::chrono::microseconds>(end_time - begin_time).count() << "[µs]" << std::endl;

    std::cout << "Solution with Gaussian elimination method is:" << std::endl;
    for (size_t i = 1; i <= dim; ++i) {
        std::cout << "x" << i << " = " << gaussian_solution.at(i, 1) << std::endl;
    }
    std::cout << "residal is:" << std::endl;
    for (size_t i = 1; i <= dim; ++i) {
        std::cout << "r" << i << " = " << gaussian_residual.at(i, 1) << std::endl;
    }

    std::cout << "Solution with Gauss-Seidel method is:" << std::endl;
    for (size_t i = 1; i <= dim; ++i) {
        std::cout << "x" << i << " = " << seidel_solution.at(i, 1) << std::endl;
    }
    std::cout << "residal is:" << std::endl;
    for (size_t i = 1; i <= dim; ++i) {
        std::cout << "r" << i << " = " << seidel_residual.at(i, 1) << std::endl;
    }

    std::cout << "Condition number of the matrix A is " << mu << std::endl;

    return 0;
}
