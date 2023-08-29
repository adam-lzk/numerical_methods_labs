#include <iostream>
#include <cmath>

using namespace std;


double f(const double& x)
{
    return sin(x) * cos(x) * cos(x);
}


double TrapezoidMethod(const double& pointA, const double& pointB, const double& accurancyEpsilon)
{
    double x = pointA;
    double trapezoidBase1;
    double trapezoidBase2;
    double integral = 0;

    while(x < pointB)
    {
        trapezoidBase1 = f(x);

        x += accurancyEpsilon;
        trapezoidBase2 = f(x);

        integral += (trapezoidBase1 + trapezoidBase2) / 2 * accurancyEpsilon;  // delta square
    }

    return integral;
}


int main()
{
    double pointA, pointB, accuracyEpsilon;
    pointA = 0;
    pointB = M_PI_2;  // M_PI_2 = pi/2
    accuracyEpsilon = pow(10, -6);  // numeric_limits<double>::min() = 2.22507e-308

    cout << TrapezoidMethod(pointA, pointB, accuracyEpsilon) << endl;
}
