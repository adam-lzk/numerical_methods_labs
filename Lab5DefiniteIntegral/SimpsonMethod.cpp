#include <iostream>
#include <cmath>

using namespace std;


double f(const double& x)
{
    return sin(x) * cos(x) * cos(x);
}


double SimpsonMethod(const double& pointA, const double& pointB, const double& accurancyEpsilon, const int& splitSegmentsQuantity)
{
    double integral = f(pointA) + f(pointB);

    int k;

    for (int i = 1; i <= splitSegmentsQuantity - 1; i++)
    {
        k = 2 + 2 * (i % 2);
        integral += k * f(pointA + i * accurancyEpsilon);

    }

    integral *= accurancyEpsilon / 3;
    return integral;
}


int main()
{
    double pointA, pointB, accuracyEpsilon;
    pointA = 0;
    pointB = M_PI_2;  // M_PI_2 = pi/2
    accuracyEpsilon = pow(10, -9);  // numeric_limits<double>::min() = 2.22507e-308

    unsigned long long splitSegmentsQuantity = 0;

    for (double i = pointA; i < pointB; i += accuracyEpsilon)
    {
        splitSegmentsQuantity++;
    }

    cout << splitSegmentsQuantity << endl;

    cout << SimpsonMethod(pointA, pointB, accuracyEpsilon, splitSegmentsQuantity) << endl;
}
