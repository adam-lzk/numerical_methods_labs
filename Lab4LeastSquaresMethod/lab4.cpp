// g++ -std=c++11 lab4.cpp -o lab4
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;


void CalculateXYsum(double& xySum, const vector<double>& givenX, const vector<double>& givenY)
{
    for(int i = 0; i < givenX.size(); i++)
    {
        xySum += givenX[i] * givenY[i];
    }
}


void CalculateElementsSum(double& elementsSum, const vector<double>& givenElements)
{
    for (const double& element : givenElements)
    {
        elementsSum += element;
    }
}


void CalculateXSquaredSum(double& xSquaredSum, const vector<double>& givenX)
{
    for (const double& element : givenX)
    {
        xSquaredSum += element * element;
    }
}


void CalculateACoef(double& aCoefficient, double xySum, double xSum, double ySum, double xSquaredSum, int vectorSize)
{
    aCoefficient = (vectorSize * xySum - xSum * ySum) / (vectorSize * xSquaredSum - xSum * xSum);
}


void CalculateBCoef(double& bCoefficient, double aCoefficient, double xSum, double ySum, int vectorSize)
{
    bCoefficient = (ySum - aCoefficient * xSum) / vectorSize;
}


void OutputInformation(const double& aCoefficient, const double& bCoefficient, const vector<double>& givenX, const vector<double>& givenY)
{
    cout << "\ngiven x:\n";
    for (const double& elements : givenX)
    {
        cout << elements << setw(8);
    }
    cout << endl << endl;

    cout << "given y:\n";
    for (const double& elements : givenY)
    {
        cout << elements << setw(8);
    }
    cout << endl << endl;

    cout << "coefficient a = " << aCoefficient << endl;
    cout << "coefficient b = " << bCoefficient << endl << endl;

    cout << "required function: y = " << aCoefficient << " x + " << bCoefficient << endl << endl;
}


int main()
{
    vector<double> givenX = {19.1, 25, 30.1, 36, 40, 45.1, 50};
    vector<double> givenY = {76.3, 77.8, 79.75, 80.8, 82.35, 83.9, 85};

    double xySum = 0;
    CalculateXYsum(xySum ,givenX, givenY);

    double xSum = 0;
    CalculateElementsSum(xSum ,givenX);

    double ySum = 0;
    CalculateElementsSum(ySum ,givenY);

    double xSquaredSum = 0;
    CalculateXSquaredSum(xSquaredSum ,givenX);

    int vectorSize = givenX.size();

    double aCoefficient;  // y = ax + b
    CalculateACoef(aCoefficient, xySum, xSum, ySum, xSquaredSum, vectorSize);

    double bCoefficient;
    CalculateBCoef(bCoefficient, aCoefficient, xSum, ySum, vectorSize);

    OutputInformation(aCoefficient, bCoefficient, givenX, givenY);

    return 0;
}
