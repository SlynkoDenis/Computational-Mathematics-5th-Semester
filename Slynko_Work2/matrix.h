#pragma once

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <ostream>
#include <set>
#include <utility>
#include <vector>

const long double pi = 2.0L * std::acos(0.0L);

// TODO: zero number of rows and columns hasn't been tested
template <typename T>
class Matrix {
private:
    std::vector<T> elements_;
    size_t rows_;
    size_t columns_;

public:
    Matrix() : rows_(0), columns_(0) {};
    Matrix(size_t tmp_rows_, size_t tmp_columns_, T initial_value = static_cast<T>(0.0L)) {
        if (tmp_rows_ == 0 || tmp_columns_ == 0)
            throw std::invalid_argument("impossible to create a matrix with 0 number of columns or rows");
        if (std::numeric_limits<size_t>::max() / tmp_rows_ <= tmp_columns_)
            throw std::invalid_argument("input dimension variables exceed size_t capacity");

        rows_ = tmp_rows_;
        columns_ = tmp_columns_;

        elements_.resize(rows_ * columns_, initial_value);
    };
    Matrix(const Matrix& other) = default;
    Matrix(Matrix&& other) noexcept = default;
    Matrix& operator=(Matrix&& other) noexcept = default;
    Matrix& operator=(const Matrix& other) = default;
    ~Matrix() noexcept = default;
    const T& at(size_t i, size_t j) const {
        if (i > rows_ || i == 0) {
            throw std::out_of_range("number of row");
        }
        if (j > columns_ || j == 0) {
            throw std::out_of_range("number of column");
        }

        return elements_[columns_ * (i - 1) + j - 1];
    };
    T& at(size_t i, size_t j) {
        if (i > rows_ || i == 0) {
            throw std::out_of_range("number of row");
        }
        if (j > columns_ || j == 0) {
            throw std::out_of_range("number of column");
        }

        return elements_[columns_ * (i - 1) + j - 1];
    };
    size_t get_rows_number() const {
        return rows_;
    };
    size_t get_columns_number() const {
        return columns_;
    };
    void swap_rows(size_t num1, size_t num2) {
        if (num1 > rows_ || num2 > rows_ || num1 == 0 || num2 == 0) {
            throw std::out_of_range("numbers of rows");
        }
        if (num1 == num2) {
            return;
        }

        std::swap_ranges(elements_.begin() + (num1 - 1) * rows_, elements_.begin() + num1 * rows_,
                         elements_.begin() + (num2 - 1) * rows_);
    }
};

namespace staff_functions {
    template <typename T>
    T naive_determinant(const Matrix<T>& obj) {
        if (obj.get_columns_number() != obj.get_rows_number())
            throw std::invalid_argument("given matrix is not square");

        size_t dim = obj.get_columns_number();
        if (dim == 2)
            return obj.at(1, 1) * obj.at(2, 2) - obj.at(1, 2) * obj.at(2, 1);

        decltype(auto) result = static_cast<T>(0.0L);
        for (size_t i = 1; i <= dim; ++i) {
            result += std::pow(-1.0, i - 1) * obj.at(1, i) * naive_determinant(form_minor_matrix(obj, 1, i));
        }

        return result;
    }

    template <typename T>
    std::pair<Matrix<T>, int> LUP_parces(const Matrix<T>& obj) {
        if (obj.get_columns_number() != obj.get_rows_number()) {
            throw std::invalid_argument("given matrix is not square");
        }

        size_t dim = obj.get_rows_number();
        Matrix<T> C(obj);
        Matrix<T> P(dim, dim);
        for (size_t i = 1; i <= dim; ++i)
            P.at(i, i) = 1.0;

        int num_of_permutations_in_P = 0;

        for (size_t i = 1; i <= dim; i++) {
            decltype(auto) pivot = static_cast<T>(0.0L);
            size_t pv_row = 0;

            for (size_t row = i; row <= dim; ++row) {
                if (std::fabs(C.at(row, i)) > pivot) {
                    pivot = std::fabs(C.at(row, i));
                    pv_row = row;
                }
            }

            if (pivot != static_cast<T>(0.0L)) {
                if (pv_row != i)
                    ++num_of_permutations_in_P;

                P.swap_rows(pv_row, i);
                C.swap_rows(pv_row, i);

                for (size_t j = i + 1; j <= dim; ++j) {
                    C.at(j, i) /= C.at(i, i);

                    const T& tmp = C.at(j, i);
                    for (size_t k = i + 1; k <= dim; ++k)
                        C.at(j, k) -= tmp * C.at(i, k);
                }
            }
            else {
                throw std::invalid_argument("the matrix is singular");
            }
        }

        return std::make_pair(C, num_of_permutations_in_P);
    }

    template <typename T>
    int is_triangular(const Matrix<T>& obj) {
        if (obj.get_rows_number() != obj.get_columns_number())
            return 0;
        size_t dim = obj.get_columns_number();

        bool flag = true;
        for (size_t i = 1; i < dim; ++i) {
            if (flag) {
                for (size_t j = i + 1; j <= dim; ++j) {
                    if (obj.at(i, j) != static_cast<T>(0.0L)) {
                        flag = false;
                        break;
                    }
                }
            }
        }
        if (flag) {
            return 1;
        }

        for (size_t i = 2; i <= dim; ++i) {
            for (size_t j = 1; j < i; ++j) {
                if (obj.at(i, j) != static_cast<T>(0.0L)) {
                    return 0;
                }
            }
        }
        return -1;
    }
}

template <typename T>
bool operator== (const Matrix<T>& lhs, const Matrix<T>& rhs) {
    if (lhs.get_rows_number() != rhs.get_rows_number() || lhs.get_columns_number() != rhs.get_columns_number())
        return false;
    for (size_t i = 1, rows_number = lhs.get_rows_number(); i <= rows_number; ++i) {
        for (size_t j = 1, columns_number = lhs.get_columns_number();
             j <= columns_number; ++j) {
            if (lhs.at(i, j) != rhs.at(i, j))
                return false;
        }
    }

    return true;
}

template <typename T>
bool operator!= (const Matrix<T>& lhs, const Matrix<T>& rhs) {
    return !(lhs == rhs);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& obj) {
    size_t rows_number = obj.get_rows_number();
    size_t columns_number = obj.get_columns_number();
    for (size_t i = 1; i < rows_number; ++i) {
        for (size_t j = 1; j <= columns_number; ++j)
            os << std::fixed << std::setprecision(12) << std::setw(3) << obj.at(i, j) << " ";
        os << std::endl;
    }
    for (size_t j = 1; j < columns_number; ++j) {
        os << std::fixed << std::setprecision(12) << std::setw(3) << obj.at(rows_number, j) << " ";
    }
    os << std::fixed << std::setprecision(12) << std::setw(3) << obj.at(rows_number, columns_number);

    return os;
}

template <typename T>
Matrix<T> operator+ (const Matrix<T>& lhs, const Matrix<T>& rhs) {
    if (lhs.get_rows_number() != rhs.get_rows_number() || lhs.get_columns_number() != rhs.get_columns_number())
        throw std::invalid_argument("size conflict of the input matrices");

    Matrix<T> result_matrix(lhs.get_rows_number(), lhs.get_columns_number());

    for (size_t i = 1, rows_number = lhs.get_rows_number(); i <= rows_number; ++i) {
        for (size_t j = 1, columns_number = lhs.get_columns_number();
             j <= columns_number; ++j) {
            result_matrix.at(i, j) = lhs.at(i, j) + rhs.at(i, j);
        }
    }

    return result_matrix;
}

template <typename T>
Matrix<T> operator- (const Matrix<T>& lhs, const Matrix<T>& rhs) {
    if (lhs.get_rows_number() != rhs.get_rows_number() || lhs.get_columns_number() != rhs.get_columns_number())
        throw std::invalid_argument("size conflict of the input matrices");

    Matrix<T> result_matrix(lhs.get_rows_number(), lhs.get_columns_number());

    for (size_t i = 1, rows_number = lhs.get_rows_number(); i <= rows_number; ++i) {
        for (size_t j = 1, columns_number = lhs.get_columns_number();
             j <= columns_number; ++j) {
            result_matrix.at(i, j) = lhs.at(i, j) - rhs.at(i, j);
        }
    }

    return result_matrix;
}

template <typename T>
Matrix<T> transpose(const Matrix<T>& other) {
    Matrix<T> result_matrix(other.get_columns_number(), other.get_rows_number());

    for (size_t i = 1, rows_number = other.get_rows_number(); i <= rows_number; ++i) {
        for (size_t j = 1, columns_number = other.get_columns_number();
             j <= columns_number; ++j) {
            result_matrix.at(j, i) = other.at(i, j);
        }
    }

    return result_matrix;
}

template <typename T>
Matrix<T> operator* (const Matrix<T>& lhs, const Matrix<T>& rhs) {
    if (lhs.get_columns_number() != rhs.get_rows_number())
        throw std::invalid_argument("size conflict of the input matrices");

    Matrix<T> result_matrix(lhs.get_rows_number(), rhs.get_columns_number());

    Matrix<T> rhs_tr = transpose(rhs);
    for (size_t i = 1, rows_number = result_matrix.get_rows_number(); i <= rows_number; ++i) {
        for (size_t j = 1, columns_number = result_matrix.get_columns_number(); j <= columns_number; ++j) {
            result_matrix.at(i, j) = static_cast<T>(0.0L);

            for (size_t k = 1; k <= rhs_tr.get_columns_number(); ++k)
                result_matrix.at(i, j) += lhs.at(i, k) * rhs_tr.at(j, k);
        }
    }

    return result_matrix;
}

template <typename T>
Matrix<T> operator* (const Matrix<T>& obj, const T& sqal) {
    Matrix<T> result_matrix(obj);

    for (size_t i = 1, rows_number = result_matrix.get_rows_number(); i <= rows_number; ++i) {
        for (size_t j = 1, columns_number = result_matrix.get_columns_number(); j <= columns_number; ++j) {
            result_matrix.at(i, j) *= sqal;
        }
    }

    return result_matrix;
}

template <typename T>
Matrix<T> operator* (const T& sqal, const Matrix<T>& obj) {
    return obj * sqal;
}

template <typename T>
Matrix<T> form_minor_matrix(const Matrix<T>& obj, size_t row_num, size_t column_num) {
    if (obj.get_rows_number() == 1 || obj.get_columns_number() == 1)
        throw std::invalid_argument("can't calculate for the inputed matrix");
    if (column_num > obj.get_columns_number() || row_num > obj.get_rows_number() ||
        column_num == 0 || row_num == 0) {
        throw std::invalid_argument("incorrect input indexes of the element");
    }

    Matrix<T> result_matrix(obj.get_rows_number() - 1, obj.get_columns_number() - 1);
    for (size_t i = 1, obj_index_row = 1, rows_number = result_matrix.get_rows_number();
         i <= rows_number; ++i, ++obj_index_row) {
        if (i == row_num) {
            ++obj_index_row;
        }

        for (size_t j = 1, obj_index_column = 1, columns_number = result_matrix.get_columns_number();
             j <= columns_number; ++j, ++obj_index_column) {
            if (j == column_num) {
                ++obj_index_column;
            }

            result_matrix.at(i, j) = obj.at(obj_index_row, obj_index_column);
        }
    }

    return result_matrix;
}


// TODO: fix
template <typename T>
Matrix<T> inverse(const Matrix<T>& obj) {
    if (obj.get_rows_number() != obj.get_columns_number())
        throw std::invalid_argument("matrix must be square");

    Matrix<T> result_matrix(obj.get_columns_number(), obj.get_columns_number());

    decltype(auto) det_obj = static_cast<T>(0.0L);

    decltype(auto) adj = static_cast<T>(0.0L);
    for (size_t i = 1, columns_number = obj.get_columns_number(); i <= columns_number; ++i) {
        for (size_t j = 1; j <= columns_number; ++j) {
            adj = std::pow(-1.0f, i + j) * determinant(form_minor_matrix(obj, j, i));
            result_matrix.at(i, j) = adj;

            if (i == 1) {
                det_obj += adj * obj.at(j, i);
            }
        }
    }

    if (det_obj == static_cast<T>(0.0L))
        throw std::invalid_argument("matrix is singular");

    result_matrix = result_matrix * (static_cast<T>(1.0L) / det_obj);

    return result_matrix;
}


// Delete this implementation
template <typename T>
Matrix<T> lower_triangular_inverse_slow(const Matrix<T>& obj) {
    if (staff_functions::is_triangular(obj) == 0) {
        throw std::invalid_argument("matrix must be triangular");
    }

    size_t dim = obj.get_columns_number();
    Matrix<T> result_matrix(dim, dim);

    for (size_t i = 1; i <= dim; ++i) {
        if (obj.at(i, i) == static_cast<T>(0.0L)) {
            throw std::invalid_argument("matrix is singular");
        }
        result_matrix.at(i, i) = static_cast<T>(1.0L) / obj.at(i, i);
    }

    decltype(auto) det = triangular_determinant(obj);
    for (size_t i = 2; i <= dim; ++i) {
        for (size_t j = 1; j < dim; ++j) {
            result_matrix.at(i, j) = std::pow(-1.0f, i + j) * determinant(form_minor_matrix(obj, j, i)) * (static_cast<T>(1.0L) / det);
        }
    }

    return result_matrix;
}

// TODO: check if the matrix is lower or higher triangular
template <typename T>
Matrix<T> lower_triangular_inverse(const Matrix<T>& obj) {
    if (staff_functions::is_triangular(obj) == 0) {
        throw std::invalid_argument("matrix must be triangular");
    }

    size_t dim = obj.get_columns_number();
    Matrix<T> result_matrix(dim, dim);

    for (size_t i = 1; i <= dim; ++i) {
        if (obj.at(i, i) == static_cast<T>(0.0L)) {
            throw std::invalid_argument("matrix is singular");
        }
        result_matrix.at(i, i) = static_cast<T>(1.0L) / obj.at(i, i);
    }

    for (size_t i = 2; i <= dim; ++i) {
        for (size_t j = i - 1; j > 0; --j) {
            T tmp = static_cast<T>(0.0L);
            for (size_t k = i; k >= j + 1; --k) {
                tmp += result_matrix.at(i, k) * obj.at(k, j);
            }

            result_matrix.at(i, j) = static_cast<T>(-1.0L) / obj.at(j, j) * tmp;
        }
    }

    return result_matrix;
}

template <typename T>
T determinant(const Matrix<T>& obj) {
    if (obj.get_rows_number() != obj.get_columns_number()) {
        throw std::invalid_argument("matrix must be square");
    }

    try {
        decltype(auto) LUP_res = staff_functions::LUP_parces(obj);
        decltype(auto) result = static_cast<T>(1.0L);
        for (size_t i = 1, rows_number = obj.get_rows_number(); i <= rows_number; ++i) {
            result *= LUP_res.first.at(i, i);
        }

        return result * static_cast<T>(std::pow(-1.0, LUP_res.second));
    } catch (std::invalid_argument& e) {
        return static_cast<T>(0.0L);
    }
}

template <typename T>
T triangular_determinant(const Matrix<T>& obj) {
    if (staff_functions::is_triangular(obj) == 0) {
        throw std::invalid_argument("matrix must be triangular");
    }

    size_t dim = obj.get_columns_number();
    T det = obj.at(1, 1);
    for (size_t i = 2; i <= dim; ++i) {
        det *= obj.at(i, i);
    }

    return det;
}

// TODO: move it to namespace
// Current implementation evaluates norm as maximum of sums of
// elements in rows taken with absolute value (so called norm 1)
template <typename T>
T norm(const Matrix<T>& A) {
    auto rows_number = A.get_rows_number();
    auto columns_number = A.get_columns_number();
    std::vector<T> sums_in_rows(rows_number, static_cast<T>(0.0L));

    for (size_t i = 1; i <= rows_number; ++i) {
        for (size_t j = 1; j <= columns_number; ++j) {
            sums_in_rows[i - 1] += std::abs(A.at(i, j));
        }
    }

    return *std::max_element(sums_in_rows.begin(), sums_in_rows.end());
}

template <typename T>
T condition_number(const Matrix<T>& A) {
    return norm(A) * norm(inverse(A));
}

// TODO: move it to namespace
template <typename T, typename UnaryPredicate, typename Condition>
std::pair<size_t, size_t> find_max_element(const Matrix<T>& A, UnaryPredicate p,
                                           Condition c) {
    size_t rows_number = A.get_rows_number();
    size_t columns_number = A.get_columns_number();

    decltype(auto) max = static_cast<T>(0.0L);
    std::pair<size_t, size_t> result = std::make_pair(0, 0);
    for (size_t j = 1; j <= columns_number; ++j) {
        if (c(j)) {
            max = p(A.at(1, j));
            result = std::make_pair(1, j);
        }
    }
    if (result.first == 0) {
        throw std::logic_error("either matrix is singular or process is finished");
    }

    for (size_t i = 1; i <= rows_number; ++i) {
        for (size_t j = 1; j <= columns_number; ++j) {
            if (c(j)) {
                decltype(auto) tmp = p(A.at(i, j));
                if (tmp > max) {
                    max = tmp;
                    result = std::make_pair(i, j);
                }
            }
        }
    }

    return result;
}

// Implemented for square matrices
template <typename T>
Matrix<T> gaussian_elimination(Matrix<T> A, Matrix<T> f) {
    size_t dim = A.get_columns_number();
    if (A.get_rows_number() != dim)
        throw std::invalid_argument("matrix A must be square");
    if (f.get_columns_number() != 1)
        throw std::invalid_argument("f must be a single column");
    if (f.get_rows_number() != dim)
        throw std::invalid_argument("f and A dimensions are different");

    std::vector<size_t> leader_indexes(dim, 0);
    std::set<size_t> indexes;
    for (size_t i = 1; i <= dim; ++i) {
        decltype(auto) leader_pos = find_max_element(A, [](const T& x){
            return std::abs(x);
        }, [indexes](size_t j){
            if (indexes.find(j) != indexes.end())
                return false;
            return true;
        });
        indexes.insert(leader_pos.second);
        leader_indexes.at(leader_pos.first - 1) = leader_pos.second;
        decltype(auto) leader = A.at(leader_pos.first, leader_pos.second);

        for (size_t j = 1; j <= dim; ++j) {
            if (j != leader_pos.first) {
                decltype(auto) coef = A.at(j, leader_pos.second) / leader;

                f.at(j, 1) -= f.at(leader_pos.first, 1) * coef;
                for (size_t k = 1; k <= dim; ++k) {
                    if (k == leader_pos.second) {
                        A.at(j, k) = static_cast<T>(0.0L);
                    } else {
                        A.at(j, k) -= A.at(leader_pos.first, k) * coef;
                    }
                }
            }
        }
    }

    decltype(auto) tmp_leader_indexes = leader_indexes;
    std::sort(tmp_leader_indexes.begin(), tmp_leader_indexes.end());
    std::vector<size_t> sequential_numbers(dim);
    std::iota(sequential_numbers.begin(), sequential_numbers.end(), 1);
    if (tmp_leader_indexes != sequential_numbers) {
        throw std::logic_error("either matrix is singular or smth went wrong with the algorithm");
    }

    Matrix<T> result(dim, 1);
    for (size_t i = 1; i <= dim; ++i) {
        result.at(leader_indexes.at(i - 1), 1) = f.at(i, 1) / A.at(i, leader_indexes.at(i - 1));
    }

    return result;
}

// Implemented for square matrices
template <typename T>
Matrix<T> gauss_seidel_method(const Matrix<T>& A, const Matrix<T>& f) {
    size_t dim = A.get_columns_number();
    if (A.get_rows_number() != dim)
        throw std::invalid_argument("matrix A must be square");
    if (f.get_columns_number() != 1)
        throw std::invalid_argument("f must be a single column");
    if (f.get_rows_number() != dim)
        throw std::invalid_argument("f and A dimensions are different");

    Matrix<T> D(dim, dim);
    for (size_t i = 1; i <= dim; ++i) {
        D.at(i, i) = A.at(i, i);
    }
    Matrix<T> L(dim, dim);
    for (size_t i = 1; i <= dim; ++i) {
        for (size_t j = 1; j < i; ++j) {
            L.at(i, j) = A.at(i, j);
        }
    }
    Matrix<T> U(dim, dim);
    for (size_t i = 1; i <= dim; ++i) {
        for (size_t j = dim; j > i; --j) {
            U.at(i, j) = A.at(i, j);
        }
    }

    decltype(auto) tmp = lower_triangular_inverse(L + D);
    decltype(auto) B = static_cast<T>(-1.0L) * tmp * U;
    decltype(auto) F = tmp * f;

    Matrix<T> iteration_result(dim, 1, static_cast<T>(1.0L));
    tmp = B * iteration_result + F;

    while (norm(iteration_result - tmp) > std::numeric_limits<T>::epsilon()) {
        iteration_result = tmp;
        tmp = B * iteration_result + F;
    }

    return iteration_result;
}

template <typename T>
std::vector<T> jacobi_eigenvalue_algorithm(Matrix<T> A) {
    if (transpose(A) != A) {
        throw std::invalid_argument("matrix must be symmetric");
    }

    size_t dim = A.get_rows_number();
    while (true) {
        T local_max = std::abs(A.at(1, 2));
        std::pair<size_t, size_t> indexes = std::make_pair(1, 2);
        for (size_t i = 1; i < dim; ++i) {
            for (size_t j = i + 1; j <= dim; ++j) {
                T tmp = std::abs(A.at(i, j));
                if (tmp > local_max) {
                    local_max = tmp;
                    indexes = std::make_pair(i, j);
                }
            }
        }

        if (local_max < std::numeric_limits<T>::epsilon()) {
            break;
        }

        long double fi = 0.25L * pi;
        if (A.at(indexes.first, indexes.first) != A.at(indexes.second, indexes.second)) {
            fi = 0.5L * std::atan(2.0L * A.at(indexes.first, indexes.second) / (A.at(indexes.first, indexes.first)
                                                                                     - A.at(indexes.second, indexes.second)));
        }
        Matrix<T> H(dim, dim);
        for (size_t i = 1; i <= dim; ++i) {
            if (i != indexes.first && i != indexes.second) {
                H.at(i, i) = static_cast<T>(1.0L);
            }
        }
        H.at(indexes.first, indexes.first) = static_cast<T>(std::cos(fi));
        H.at(indexes.second, indexes.second) = static_cast<T>(std::cos(fi));
        H.at(indexes.second, indexes.first) = static_cast<T>(std::sin(fi));
        H.at(indexes.first, indexes.second) = static_cast<T>(-1.0L * std::sin(fi));

        A = transpose(H) * A * H;
    }

    std::vector<T> result(dim, static_cast<T>(0.0L));
    for (size_t i = 1; i <= dim; ++i) {
        result.at(i - 1) = A.at(i, i);
    }

    return result;
}

template <typename T>
T operator+ (const T& lhs, const Matrix<T>& rhs) {
    if (rhs.get_columns_number() != 1 || rhs.get_rows_number() != 1)
        throw std::logic_error("can't operate");

    return lhs + rhs.at(1, 1);
}

template <typename T>
T operator- (const T& lhs, const Matrix<T>& rhs) {
    if (rhs.get_columns_number() != 1 || rhs.get_rows_number() != 1)
        throw std::logic_error("can't operate");

    return lhs - rhs.at(1, 1);
}

template <typename T>
T operator- (const Matrix<T>& lhs, const T& rhs) {
    if (lhs.get_columns_number() != 1 || lhs.get_rows_number() != 1)
        throw std::logic_error("can't operate");

    return lhs.at(1, 1) - rhs;

}
