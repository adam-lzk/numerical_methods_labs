//#include "stdafx.h" (not my code) g++ -std=c++11 lab4v4.cpp -o lab4v4
//#include <tchar.h>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;
int main()
{
    setlocale(LC_ALL, "rus");
    double* a = NULL, * b = NULL, ** sum = NULL;
    const int N = 8;

    double  x[N] = { 6.1, 7.5, 8.88, 11.1, 12.2 },
        y[N] = { 3185, 5390, 6860, 10045, 12740};
    int K, i, j, k, m;
    double z, c;
    cout << "Порядок k= ";
    cin >> K;
    cin.get();
    b = new double[K + 1];
    a = new double[K + 1];
    sum = new double* [K + 1];
    for (int i = 0; i < K + 1; i++)
        sum[i] = new double[K + 1];



    for (i = 0; i < K + 1; i++)
    {
        for (j = 0; j < K + 1; j++)
        {
            sum[i][j] = 0;
            for (k = 0; k < N; k++)
            {
                sum[i][j] += pow(x[k], i + j);
            }
        }
    }
    for (i = 0; i < K + 1; i++)
    {
        b[i] = 0;
        for (k = 0; k < N; k++)
        {
            b[i] += pow(x[k], i) * y[k];
        }
    }
    for (int i = 0; i < K + 1; i++)
    {
        for (int j = 0; j < K + 1; j++)
            cout << sum[i][j] << "\t";
        cout << endl;
    }

    for (int i = 0; i < K + 1; i++)
        cout << b[i] << "\t";
    cout << endl;
    for (int i = 0; i < K; i++)
    {
        m = i;
        for (int j = i + 1; j < K + 1; j++)
        {
            if (fabs(sum[m][i]) < fabs(sum[j][i]))
                m = j;
        }
        for (int k = i; k < K + 1; k++)
        {
            z = sum[m][k];
            sum[m][k] = sum[i][k];
            sum[i][k] = z;
        }
        z = b[m];  b[m] = b[i];  b[i] = z;
        for (int i = 0; i < K + 1; i++)
        {
            for (int j = 0; j < K + 1; j++)
            {
                cout << sum[i][j] << setw(15);
            }
            cout << "b" << i << "= ";
            cout << b[i] << endl;
        }
        cin.get();

        for (int j = i + 1; j < K + 1; j++)
        {
            c = -sum[j][i] / sum[i][i];
            cout << "\n\n ! " << c << endl;
            for (int k = i; k < K + 1; k++)
            {
                sum[j][k] = sum[j][k] + c * sum[i][k];
            }
            b[j] = b[j] + c * b[i];
        }
    }

    for (int i = 0; i < K + 1; i++)
    {
        for (int j = 0; j < K + 1; j++)
        {
            cout << sum[i][j] << setw(15);
        }
        cout << "b" << i << "=";
        cout << b[i] << endl;
    }
    cin.get();

    a[K] = b[K] / sum[K][K];

    for (int i = K + 1 - 2; i >= 0; i--)
    {
        for (int k = i + 1; k < K + 1; k++)
        {
            b[i] = b[i] - a[k] * sum[i][k];
        }
        a[i] = b[i] / sum[i][i];
    }


    for (int i = 0; i < K + 1; i++)
        cout << "X" << i << "=" << a[i] << endl;
    cin.get();

    return 0;
}
