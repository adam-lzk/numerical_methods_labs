// g++ -std=c++11 lab3v2.cpp -o lab3v2 (not my code)

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <algorithm>

using namespace std;

const size_t outwidth = 20;
const size_t precision = 10;
const double round_error = 0.000000001;

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

// second task

double first_equasion_1(const vector<double>& x, const int& a, const int& k) {
    return ((k - a) * x.at(1) * x.at(2)) / a;
}

double second_equasion_1(const vector<double>& x, const int& a, const int& k) {
    return ((k + a) * x.at(0) * x.at(2)) / k;
}

double third_equasion_1(const vector<double>& x, const int& a, const int& k) {
    return ((a - k) * x.at(0) * x.at(1)) / a;
}

double first_equasion(const vector<double>& x, const vector<double>& lambda) {
    return ((2 * lambda.at(0) + 4 * lambda.at(1)) * x.at(0) + (2 * lambda.at(0) - 2 * lambda.at(1)) * x.at(1) + (2 * lambda.at(0) - 2 * lambda.at(1)) * x.at(2)) / 6 + (4 * lambda.at(0) + 2 * lambda.at(1)) / 6;
}

double second_equasion(const vector<double>& x, const vector<double>& lambda) {
    return ((2 * lambda.at(0) - 2 * lambda.at(1)) * x.at(0) + (2 * lambda.at(0) + lambda.at(1) + 3 * lambda.at(2)) * x.at(1) + (2 * lambda.at(0) + lambda.at(1) - 3 * lambda.at(2)) * x.at(2)) / 6 + (4 * lambda.at(0) - lambda.at(1) - 9 * lambda.at(2)) / 6;
}

double third_equasion(const vector<double>& x, const vector<double>& lambda) {
    return ((2 * lambda.at(0) - 2 * lambda.at(1)) * x.at(0) + (2 * lambda.at(0) + lambda.at(1) - 3 * lambda.at(2)) * x.at(1) + (2 * lambda.at(0) + lambda.at(1) + 3 * lambda.at(2)) * x.at(2)) / 6 + (4 * lambda.at(0) - lambda.at(1) + 9 * lambda.at(2)) / 6;
}

double diff_first_x1(const vector<double>& lambda) {
    return (2 * lambda.at(0) + 4 * lambda.at(1)) / 6;
}

double diff_first_x2(const vector<double>& lambda) {
    return (2 * lambda.at(0) - 2 * lambda.at(1)) / 6;
}

double diff_first_x3(const vector<double>& lambda) {
    return (2 * lambda.at(0) - 2 * lambda.at(1)) / 6;
}

double diff_second_x1(const vector<double>& lambda) {
    return (2 * lambda.at(0) - 2 * lambda.at(1)) / 6;
}

double diff_second_x2(const vector<double>& lambda) {
    return (2 * lambda.at(0) + lambda.at(1) + 3 * lambda.at(2)) / 6;
}

double diff_second_x3(const vector<double>& lambda) {
    return (2 * lambda.at(0) + lambda.at(1) - 3 * lambda.at(2)) / 6;
}

double diff_third_x1(const vector<double>& lambda) {
    return (2 * lambda.at(0) - 2 * lambda.at(1)) / 6;
}

double diff_third_x2(const vector<double>& lambda) {
    return (2 * lambda.at(0) + lambda.at(1) - 3 * lambda.at(2)) / 6;
}

double diff_third_x3(const vector<double>& lambda) {
    return (2 * lambda.at(0) + lambda.at(1) + 3 * lambda.at(2)) / 6;
}

vector<vector<double>> find_jacobian_1(const vector<double>& lambda)
{
    vector<vector<double>> jacobian_1(lambda.size(), vector<double>(lambda.size(), 0));
    jacobian_1.at(0).at(0) = diff_first_x1(lambda);
    jacobian_1.at(0).at(1) = diff_first_x2(lambda);
    jacobian_1.at(0).at(2) = diff_first_x3(lambda);
    jacobian_1.at(1).at(0) = diff_second_x1(lambda);
    jacobian_1.at(1).at(1) = diff_second_x2(lambda);
    jacobian_1.at(1).at(2) = diff_second_x3(lambda);
    jacobian_1.at(2).at(0) = diff_third_x1(lambda);
    jacobian_1.at(2).at(1) = diff_third_x2(lambda);
    jacobian_1.at(2).at(2) = diff_third_x3(lambda);
    return jacobian_1;
}

// dx1
double find_d1(vector<double> x, vector<double> y, const vector<double>& lambdas)
{
    vector<double> elements = { abs(first_equasion(x, lambdas)),
                                     abs(second_equasion(x, lambdas)),
                                     abs(third_equasion(x, lambdas)) };

    return *max_element(elements.begin(), elements.end());
}

double find_d2(vector<double> x, vector<double> x1)
{
    size_t n = x.size();
    double max = 0;
    for (size_t i = 0; i < n; ++i) {
        if (abs(x1.at(i)) < 1) {
            if (abs(x1.at(i) - x.at(i)) > max) {
                max = abs(x1.at(i) - x.at(i));
            }
        }
        else {
            if (abs((x1.at(i) - x.at(i)) / x1.at(i)) > max) {
                max = abs((x1.at(i) - x.at(i)) / x1.at(i));
            }
        }
    }

    return max;
}

vector<double> solve_system_with_Newtons_method_three_equations(const double& e, vector<double>& x, const vector<double>& y, const vector<double>& lambdas, const double& h, const size_t& number_of_itterations = 500)
{
    size_t n = x.size();
    size_t k_ = 1;
    double d1 = 1, d2 = 1;
    vector<double> F(3, 0);
    vector<vector<double>> J;
    vector<double> deltax;
    vector<double> x1(3, 0);
    while ((d1 > e || d2 > e) && k_ < number_of_itterations) {
        //find_residual_vector
        F.at(0) = -(first_equasion(x, lambdas));
        F.at(1) = -(second_equasion(x, lambdas));
        F.at(2) = -(third_equasion(x, lambdas));
        //find_jacobian
        J = find_jacobian_1(lambdas);
        vector<vector<double>> rectangle = fill_rectangle_matrix(J, F);
        deltax = solve_system_with_gauss_method(rectangle);
        for (size_t i = 0; i < n; ++i) {
            x1.at(i) = x.at(i) + deltax.at(i);
        }
        d1 = find_d1(x, y, lambdas);
        d2 = find_d2(x, x1);
        //cout << setw(outwidth) << setprecision(precision) << k_ << "   " << d1 << "   " << d2 << endl;
        k_++;
        x = x1;
        if (k_ >= number_of_itterations) {
            cout << "ERROR: IER = 2" << endl;
        }
    }
    return x;
}

double find_step(double e, double hmax, vector<double>& answer) {
    double h1, h2, h3, hmin;
    h1 = e / (abs(answer.at(0)) + (e / hmax));
    h2 = e / (abs(answer.at(1)) + (e / hmax));
    h3 = e / (abs(answer.at(2)) + (e / hmax));
    if (h1 < h2) {
        hmin = h1;
    }
    else {
        hmin = h2;
    }
    if (h3 < hmin) {
        hmin = h3;
    }
    return hmin;
}

vector<double> solve_system_of_three_differential_equations_with_Eulers_explicit_method(const int& a = 1, const int& k = 2, const double& e = pow(10, -2)) {
    double t = 0, T = 1;
    vector<double> answer(3, 1);
    double hmax = 1, h;
    vector<double> answer_1(3, 0);
    size_t i = 0;

    do {
        i++;
        // initialize new answer vector with values of answer put in system of differential functions
        answer_1.at(0) = first_equasion_1(answer, a, k);
        answer_1.at(1) = second_equasion_1(answer, a, k);
        answer_1.at(2) = third_equasion_1(answer, a, k);
        // find step
        h = find_step(e, hmax, answer_1);
        // y += tf
        answer.at(0) += h * answer_1.at(0);
        answer.at(1) += h * answer_1.at(1);
        answer.at(2) += h * answer_1.at(2);
        t += h;
        cout << setprecision(6) << setw(10) /*<< i << ";" << "h = " << h << ";"*/ << t << ";"  << answer.at(0) << ";" << answer.at(1) << ";" << answer.at(2) << endl;
    } while (t < T);
    return answer;
}

double find_step_1(double h, double e, const vector<double>& ek, bool v)
{
    vector<double> hh(3, 0);
    if (v) {
        for (size_t i = 0; i < 3; ++i) {
            hh.at(i) = pow(e / abs(ek.at(i)), 0.5) * h;
        }
    }
    else {
        for (size_t i = 0; i < 3; ++i) {
            if (abs(ek.at(i)) > e) {
                hh.at(i) = h / 2;
            }
            if (e / 4 < abs(ek.at(i)) && abs(ek.at(i)) <= e) {
                hh.at(i) = h;
            }
            if (abs(ek.at(i)) < e / 4) {
                hh.at(i) = 2 * h;
            }
        }
    }
    return *min_element(hh.begin(), hh.end());
}


vector<double> find_new_eps(vector<double> uk_1, vector<double> u, vector<double> uk1, double h, double hk_1)
{
    vector<double> ek(3, 0);
    ek.at(0) = -(h / (h + hk_1)) * (uk1.at(0) - u.at(0) - (h * (u.at(0) - uk_1.at(0))) / (hk_1));
    ek.at(1) = -(h / (h + hk_1)) * (uk1.at(1) - u.at(1) - (h * (u.at(1) - uk_1.at(1))) / (hk_1));
    ek.at(2) = -(h / (h + hk_1)) * (uk1.at(2) - u.at(2) - (h * (u.at(2) - uk_1.at(2))) / (hk_1));
    return ek;
}

vector<double> solve_system_of_three_differential_equations_with_Eulers_implicit_method(const vector<double>& start_x, const vector<double>& lambdas, const double& T, const bool& mode = 0, const double& e = pow(10, -5)) {
    double t = 0;
    vector<double> answer = start_x;
    double hmax = 1;
    double hmin = 0.01;
    double h = hmin;
    double hk1 = hmin;
    double hk_1 = hmin;
    double tk;
    vector<double> uu(3, 0);
    // vectors for saving answer in another variable
    vector<double> uk1(3, 0);
    vector<double> uk_1(3, 0);

    vector<double> ek(3, 0);
    int i = 0, z;
    do {
        do {
            z = 0;
            i++;
            tk = t + h;
            // solve system with newton's method
            uk1 = solve_system_with_Newtons_method_three_equations(e, uk1, answer, lambdas, h);
            // find eps with formula
            ek = find_new_eps(uk_1, answer, uk1, h, hk_1);
            if (abs(ek.at(0)) > e || abs(ek.at(1)) > e || abs(ek.at(2)) > e) {
                h /= 2;
                hk1 = h;
                uk1.at(0) = answer.at(0);
                uk1.at(1) = answer.at(1);
                uk1.at(2) = answer.at(2);
                z = 1;
            }
        } while (z == 1);
        // the last coefficient defines one of two methods: 1 - quasi-optimal, 0 - method of three zones
        hk1 = find_step_1(h, e, ek, mode);
        if (hk1 > hmax) {
            hk1 = hmax;
        }
        uk_1.at(0) = answer.at(0);
        uk_1.at(1) = answer.at(1);
        uk_1.at(2) = answer.at(2);
        answer.at(0) = uk1.at(0);
        answer.at(1) = uk1.at(1);
        answer.at(2) = uk1.at(2);
        hk_1 = h;
        h = hk1;
        t = tk;
        cout << setprecision(6) << fixed << t << ";" << answer.at(0) << ";" << answer.at(1) << ";" << answer.at(2) << endl;
    } while (t < T);
    return answer;
}

int main()
{
    vector<double> start_x = { 10, 22, 9 };
    vector<double> lambdas = { -1, -1, -1 };
    print_vector(solve_system_of_three_differential_equations_with_Eulers_explicit_method());
    cout << "\n\n\n\n";
    print_vector(solve_system_of_three_differential_equations_with_Eulers_implicit_method(start_x, lambdas, 1));
    return 0;
}
