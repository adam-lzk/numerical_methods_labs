// компиляция - g++ -std=c++11 lab4v2.cpp -o lab4v2 (not my code)

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

const size_t outwidth = 20;
const double round_error = pow(10, -6);

template<typename T>
void print_matrix(const vector<vector<T>>& matrix) {
    const size_t rows_num = matrix.size();
    const size_t column_num = matrix.at(0).size();
    for (size_t i = 0; i < rows_num; ++i) {
        for (size_t j = 0; j < column_num; ++j) {
            cout << setw(outwidth) << matrix.at(i).at(j);
        }
        cout << endl;
    }
}

template<typename T>
void print_vector(const vector<T>& const_vector) {
    size_t length = const_vector.size();
    for (size_t i = 0; i < length; ++i) {
        cout << setw(outwidth) << const_vector.at(i);
    }
    cout << endl;
}

template<typename T>
vector<vector<T>> fill_rectangle_matrix(const vector<vector<T>>& matrix, const vector<T>& const_vector) {
    vector<vector<T>> solvation_matrix;
    const size_t rows_num = matrix.size();
    // filling matrix with square matrix and vector as last column
    for (size_t i = 0; i < rows_num; ++i) {
        vector<T> row;
        for (size_t j = 0; j <= rows_num; ++j) {
            if (j < rows_num) {
                row.push_back(matrix.at(i).at(j));
            }
            else {
                row.push_back(const_vector.at(i));
            }
        }
        solvation_matrix.push_back(row);
    }
    return solvation_matrix;
}

template<typename T>
void swap_rows_for_not_null_diagonal(vector<vector<T>>& matrix) {
    size_t matrix_size = matrix.size();
    for (size_t i = 0; i < matrix_size; ++i) {
        if (abs(matrix.at(i).at(i)) <= round_error) {
            if (i > 0) {
                swap(matrix.at(i), matrix.at(i - 1));
            }
            else {
                swap(matrix.at(i), matrix.at(i + 1));
            }
        }
    }
}

template<typename T>
vector<T> solve_system_with_gauss_method(vector<vector<T>>& solvation, const int& epsilon = 0) {
    swap_rows_for_not_null_diagonal(solvation);
    const size_t heigth = solvation.size();
    const size_t width = solvation.at(0).size();
    // forward pass in gauss method
    for (size_t i = 0; i < heigth; ++i) {
        // get diagonal element
        T aii = solvation.at(i).at(i);
        // divide row by it's first element (not zero in this algorithm)
        for (size_t j = i; j < width; ++j) {
            solvation.at(i).at(j) /= aii;
        }
        // perform step subtract one equation to next
        for (size_t j = i + 1; j < heigth; ++j) {
            T aji = solvation.at(j).at(i);
            for (size_t k = i + 1; k < width; ++k) {
                solvation.at(j).at(k) -= solvation.at(i).at(k) * aji;
            }
            // made zerolike triangle bottom in matrix
            solvation.at(j).at(i) = 0;
        }
    }
    // backward pass in gauss method
    for (int i = width - 2; i > 0; --i) {
        // perform step subtract one equation to next
        for (int j = i - 1; j >= 0; --j) {
            T aji = solvation.at(j).at(i);
            for (int k = width - 1; k >= i; --k) {
                solvation.at(j).at(k) -= solvation.at(i).at(k) * aji;
            }
        }
    }
    vector<T> answer;
    if (epsilon) {
        for (size_t i = 0; i < heigth; ++i) {
            answer.push_back((round((solvation.at(i).at(width - 1)) * pow(10, epsilon))) / pow(10, epsilon));
        }
    }
    else {
        for (size_t i = 0; i < heigth; ++i) {
            answer.push_back(solvation.at(i).at(width - 1));
        }
    }
    return answer;
}


// N - number of experiments, m - power of approximating polynomial
template<typename T>
vector<T> minimal_square_method(const vector<T>& x, vector<T>& y, const size_t& M) {
    size_t m = M + 1;
    size_t N = x.size();
    vector<T> POWERX;
    for (size_t i = 1; i <= 2 * (m-1); ++i) {
        T staff = 0;
        for (size_t j = 0; j < N; ++j) {
            staff += pow(x.at(j), i);//N-ный икс возводим в i-тую степень и суммируем
        }
        POWERX.push_back(staff);
    }
    cout << "POWERX:\n";
    print_vector(POWERX);
    vector<vector<T>> SUMX(m, vector<T>(m, 0));
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < m; ++j) {
            if (!i && !j) {
                SUMX.at(0).at(0) = N;
            }
            else {
                SUMX.at(i).at(j) = POWERX.at(i + j - 1);
            }
        }
    }
    cout << "\nSUMX:\n";
    print_matrix(SUMX);
    cout << "\n";
    vector<T> right_part;
    for (size_t i = 0; i < m; ++i) {
        T staff = 0;
        for (size_t j = 0; j < N; ++j) {
            staff += y.at(j) * pow(x.at(j), i);
        }
        right_part.push_back(staff);
    }
    cout << "\nright_part:\n";
    print_vector(right_part);
    vector<vector<T>> rectangle = fill_rectangle_matrix(SUMX, right_part);
    vector<T> answers = solve_system_with_gauss_method(rectangle);
    T remaining_variance = 0;
    for (size_t i = 0; i < N; ++i) {//сумма от 1 до N (yi - a0 -a1xi)
        T staff = y.at(i);
        for (size_t j = 0; j < m; ++j) {
            staff -= answers.at(j) * pow(x.at(i), j);
        }
        remaining_variance += pow(staff, 2);
    }
    cout << " \nremaining variance = " << remaining_variance / (N - m - 1) << endl; //остаточная дисперсия
    cout << " \nstandard deviation = " << sqrt(remaining_variance / (N - m - 1)) << endl; //стандартное отклонение
    return answers;
}

int main()
{
    vector<double> x = { 1.1, 1.4, 1.7, 2.1, 2.6, 4.7, 6.1, 7.0, 10.0, 12.8, 16.5, 20.8, 40.6 };
    vector<double> y = { 25.0, 22.7, 22.1, 19.8, 17.0, 12.3, 10.7, 10.0, 8.2, 6.7, 5.6, 5.0, 3.5 };
    size_t n = x.size();
    for (size_t i = 0; i < n; ++i) {
        y[i] = log(y[i]);
        x[i] = log(x[i]);
    }
    vector<double> coef = minimal_square_method(x, y, 1);
    coef[1] = -1 / coef[1];
    coef[0] = exp(coef[0]);
    print_vector(coef);
    system("pause");
    return 0;
}