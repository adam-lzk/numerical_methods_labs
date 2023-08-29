// Ax = B slough solution by the Gauss method (not my algorithm)

#include <iostream>
#include <iomanip>

using namespace std;


double** CreateArray(int order)
{
    double** matrix = new double* [order];

    for (int i = 1; i <= order; i++)
        matrix[i] = new double [order];

    return matrix;
}


void FillArray(double** matrixA, double** matrixAcopy, double* vectorB, double *vectorBcopy, int order)
{
    for (int i = 1; i <= order; i++)
    {
        for (int j = 1; j <= order; j++)
        {
            cout << "a[" << i << "][" << j << "] = ";
            cin >> matrixA[i][j];
            matrixAcopy[i][j] = matrixA[i][j];
        }

        cout << "b[" << i << "] = ";
        cin >> vectorB[i];
        vectorBcopy[i] = vectorB[i];
    }
}


bool isLinearlyDependent(double** matrixA, int order)
{
    // ...
    return true;
}


void SelectingReferenceElement(double** matrixA, double* vectorB, int order)
{
    double* tempVector = new double [order];
    double tempElement;

    for (int n = 1; n <= order; n++)
    {
        for (int j = n + 1; j <= order; j++)
        {
            if (matrixA[n][n] == 0)  // проверка опорного элемента на == 0, если == 0, то ищем в этом же столбце ненулевой элемент и свапаем эти строки
            {
                for (int i = n + 1; i <= order; i++)
                {
                    if (matrixA[i][n] != 0)
                    {
                        for (int p = 1; p <= order; p++)  // swap strings
                        {
                            tempVector[p] = matrixA[n][p];
                            matrixA[n][p] = matrixA[i][p];
                            matrixA[i][p] = tempVector[p];

                            tempElement = vectorB[n];
                            vectorB[n] = vectorB[i];
                            vectorB[i] = tempElement;
                        }
                        delete [] tempVector;
                        break;
                    }
                }
            }
        }
    }
}


void ForwardSubstitution(double** matrixA, double* vectorB, double* vectorX, int order, double temp1, double temp2)
{
    for (int n = 1; n <= order; n++)  // прямой ход
    {
        for (int j = n + 1; j <= order; j++)
        {
            SelectingReferenceElement(matrixA, vectorB, order);

            temp1 = matrixA[j][n] / matrixA[n][n];  // формула (1)

            for (int i = n; i <= order; i++)
            {
                matrixA[j][i] -= temp1 * matrixA[n][i];  // формула (2)
            }

            vectorB[j] -= temp1 * vectorB[n];  // формула (3)
        }
    }
}


void BackSubstitution(double** matrixA, double* vectorB, double* vectorX, int order, double temp1, double temp2)
{
    for (int n = order; n >= 1; n--)  // обратный ход
    {
        temp1 = 0;

        for (int j = n + 1; j <= order; j++)
        {
            temp2 = matrixA[n][j] * vectorX[j];  // формула (4)
            temp1 = temp1 + temp2;  // формула (4)
        }

        vectorX[n] = (vectorB[n] - temp1) / matrixA[n][n];  // формула (4)
    }
}


void GaussMethodSolution(double** matrixA, double* vectorB, double* vectorX, int order)  // выбор максимального элемента!!!!!!!!!!!!!!!!!!!!!!
{
    double temp1, temp2;

    ForwardSubstitution(matrixA, vectorB, vectorX, order, temp1, temp2);

    BackSubstitution(matrixA, vectorB, vectorX, order, temp1, temp2);
}


void DiscrepancyVectorCalculation(double** matrixAcopy, double* vectorB, double* vectorBcopy, double* vectorX, double *discrepancyVector, int order)
{
    for (int i = 1; i <= order; i++)
    {
        double sum = 0;

        for (int j = 1; j <= order; j++)
        {
            sum += matrixAcopy[i][j] * vectorX[j];
        }

        discrepancyVector[i] = sum - vectorBcopy[i];
    }
}


void OutputArray(double** matrixA, double* vectorB, int order)
{
    cout << endl;

    for (int i = 1; i <= order; i++)
    {
        for (int j = 1; j <= order; j++)
        {
            cout << setw(10) << matrixA[i][j];
        }
        cout << setw(6) << "|" << setw(10) << vectorB[i] << endl << endl << endl;
    }
}


void DeleteAllArrays(double** matrixA, double** matrixAcopy, double* vectorB, double* vectorBcopy, double* vectorX, double *discrepancyVector, int order)
{
    for (int i = 1; i <= order; i++)
        delete [] matrixA[i];

    delete [] matrixA;


    for (int i = 1; i <= order; i++)
        delete [] matrixAcopy[i];

    delete [] matrixAcopy;


    delete [] vectorB;

    delete [] vectorX;

    delete [] discrepancyVector;
}


int main()
{
    int order;
    cout << "order = ";
    cin >> order;

    double** matrixA = CreateArray(order);
    double** matrixAcopy = CreateArray(order);

    double *vectorB = new double[order];
    double *vectorBcopy = new double[order];

    double *vectorX = new double[order];

    FillArray(matrixA, matrixAcopy, vectorB, vectorBcopy, order);

    OutputArray(matrixA, vectorB, order);

    GaussMethodSolution(matrixA, vectorB, vectorX, order);

    cout << "  the roots of the equation:\n";
    for (int i = 1; i <= order; i++)
    {
        cout << "x[" << i << "] = " << vectorX[i] << endl;
    }


    double *discrepancyVector = new double[order];

    DiscrepancyVectorCalculation(matrixAcopy, vectorB, vectorBcopy, vectorX, discrepancyVector, order);

    cout << "\n\n  the discrepancy vector:\n";
    // cout << fixed << setprecision(8);
    for (int i = 1; i <= order; i++)
    {
        cout << "dv[" << i <<"] = " << discrepancyVector[i] << endl;  // ???
    }
    cout << fixed << setprecision(4);

    cout << "\n\n  upper triangular matrix A:\n";
    OutputArray(matrixA, vectorB, order);

    DeleteAllArrays(matrixA, matrixAcopy, vectorB, vectorBcopy, vectorX, discrepancyVector, order);

    return 0;
}
