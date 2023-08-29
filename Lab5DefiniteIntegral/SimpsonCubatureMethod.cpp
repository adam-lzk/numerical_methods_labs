// компиляция - g++ -std=c++11 SimpsonCubatureMethod.cpp -o SimpsonCubatureMethod (not my algorithm)

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


double f(const double& x, const double& y)
{
    return (pow(x, 2) + 2 * y);
}


double SimpsonCubatureMethod(const vector<double>& first_span, const vector<double>& second_span, const double& accuracy)
{
    int n = 4, m = n;
    double step1 = (first_span.at(1) - first_span.at(0)) / (2 * n), step2 = (second_span.at(1) - second_span.at(0)) / (2 * m);
    double integral = 0, integral2;

    vector<double> x = { first_span.at(0) };
    vector<double> y = { second_span.at(0) };

    for (int i = 1; i <= 2 * (n + 2); ++i) {
        x.push_back(first_span.at(0) + i * step1);
        y.push_back(second_span.at(0) + i * step2);
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            integral += f(x.at(2 * i), y.at(2 * j)) + f(x.at(2 * i + 2), y.at(2 * j)) + f(x.at(2 * i + 2), y.at(2 * j + 2)) + f(x.at(2 * i), y.at(2 * j + 2)) +
                4 * (f(x.at(2 * i + 1), y.at(2 * j)) + f(x.at(2 * i + 2), y.at(2 * j + 1)) + f(x.at(2 * i + 1), y.at(2 * j + 2)) + f(x.at(2 * i), y.at(2 * j + 1))) +
                16 * f(x.at(2 * i + 1), y.at(2 * j + 1));
        }
    }
    integral *= step1 * step2 / 9;

    do {
        n *= 2;
        step1 = (first_span.at(1) - first_span.at(0)) / (2 * n);
        step2 = (second_span.at(1) - second_span.at(0)) / (2 * m);
        integral2 = integral;
        integral = 0;
        vector<double> x = { first_span.at(0) };
        vector<double> y = { second_span.at(0) };

        for (int i = 1; i <= 2 * (n + 2); ++i) {
            x.push_back(first_span.at(0) + i * step1);
            y.push_back(second_span.at(0) + i * step2);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                integral += f(x.at(2 * i), y.at(2 * j)) + f(x.at(2 * i + 2), y.at(2 * j)) + f(x.at(2 * i + 2), y.at(2 * j + 2)) + f(x.at(2 * i), y.at(2 * j + 2)) +
                    4 * (f(x.at(2 * i + 1), y.at(2 * j)) + f(x.at(2 * i + 2), y.at(2 * j + 1)) + f(x.at(2 * i + 1), y.at(2 * j + 2)) + f(x.at(2 * i), y.at(2 * j + 1))) +
                    16 * f(x.at(2 * i + 1), y.at(2 * j + 1));
            }
        }
        integral *= step1 * step2 / 9;
    } while (abs(integral - integral2) > 15 * accuracy);

    cout << integral;
    return integral;
}


int main()
{
    double first_span = 3;
    double second_span = 4.254; 
    vector<double> first_span_1 = { 0, 2.0 };
    vector<double> second_span_1 = { 0, 1.0 };
    double epsilon_1 = pow(10, -4);
    double epsilon_2 = pow(10, -5);

    cout << "\nSimpson's cubature method:\n\n";
    cout << "epsilon 10^(-4):\n";
    SimpsonCubatureMethod(first_span_1, second_span_1, epsilon_1);

    cout << "\n\nepsilon 10^(-5):\n";
    SimpsonCubatureMethod(first_span_1, second_span_1, epsilon_2);

    cout << endl << endl;
    return 0;
}
