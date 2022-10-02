// Ax = B slough solution by the Jordan-Gauss method
// hi

#include <iostream>
#include <iomanip>

using namespace std;


double** CreateArray(int order)
{
    double** matrixA = new double* [order];

    for (int i = 0; i < order; i++)
    {
        matrixA[i] = new double [order];
    }

    return matrixA;
}


void FillArray(double** matrixA, double* vectorB, int order)
{
    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {
            cout << "a[" << i+1 << "][" << j+1 << "] = ";
            cin >> matrixA[i][j];
        }

        cout << "b[" << i+1 << "] = ";
        cin >> vectorB[i];
    }
}


double* JordanGaussMethod(double** matrixA, double* vectorB, int order)
{
    for (int repeat = 0; repeat < order; repeat++)
    {
        for (int i = 0; i < order; i++)  // i = ? (repeat + 1)
        {
            // if matrixA[i][i] != 0
            for (int j = repeat + 1; j < order; j++)  // делю i-ю строку расширенной матрицы на значение элемента matrixA[i][repeat] (первого ненулевого элемента)
            {
                matrixA[i][j] /= matrixA[i][repeat];
            }
            vectorB[i] /= matrixA[i][repeat];

            matrixA[i][repeat] = 1;
        }

        for (int i = repeat + 1; i < order; i++)
        {
            for (int j = repeat; j < order; j++)
            {
                matrixA[i][j] -= matrixA[repeat][j];
            }
            vectorB[i] -= vectorB[repeat];
        }
    }

    vectorB[order - 1] /= matrixA[order - 1][order - 1];

    matrixA[order - 1][order - 1] = 1;

    // the matrix has been reduced to a upper triangular form

    ArrayOutput(matrixA, vectorB, order);

    double *vectorX = new double[order];

    for (int reverseRpeat = order - 1; reverseRpeat >= 0; reverseRpeat--)  // reverseRpeat (order - 1  ->  0)
    {
        vectorX[reverseRpeat] = vectorB[reverseRpeat];

        for (int i = reverseRpeat - 1; i >= 0; i--)
        {
            for (int j = reverseRpeat; j >= 0; j--)
            {
                vectorX[reverseRpeat] = vectorB[reverseRpeat] - matrixA[i][j] * vectorX[reverseRpeat];
            }
        }

        


    }

    return vectorX;
}


void ArrayOutput(double** matrixA, double* vectorB, int order)
{
    cout << endl;

    for (int i = 0; i < order; i++)
    {
        for (int j = 0; j < order; j++)
        {
            cout << setw(5) << matrixA[i][j];
        }
        cout << setw(3) << "|" << setw(3) << vectorB[i] << endl << endl;
    }
}


void DeleteArray(double** matrixA, double* vectorB, int order)
{
    for (int i = 0; i < order; i++)
        delete [] matrixA[i];

    delete [] matrixA;

    delete [] vectorB;
}


int main()
{
    int order;
    cout << "order = ";
    cin >> order;

    double** matrixA = CreateArray(order);
    double *vectorB = new double[order];

    FillArray(matrixA, vectorB, order);

    JordanGaussMethod(matrixA, vectorB, order);

    ArrayOutput(matrixA, vectorB, order);

    DeleteArray(matrixA, vectorB, order);

    return 0;
}
