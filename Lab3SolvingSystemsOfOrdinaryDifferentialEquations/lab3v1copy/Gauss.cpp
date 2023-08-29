#pragma once
#include "Gauss.h"

double** createMatrix(int n)
{
	double **matrix = new double*[n];
	for (int i = 0; i < n; i++)
		matrix[i] = new double[n];
	return matrix;
}
void deleteMatrix(double **matrix, int n)
{
	for (int i = 0; i < n; i++)
		delete matrix[i];
	delete[] matrix;
}
double* createVector(int n)
{
	double *vector = new double[n];
	return vector;
}
void deleteVector(double* vector)
{
	delete[] vector;
}
void copy(double* old, double* newvec, const int n)
{
	for (int i = 0; i < n; ++i)
		newvec[i] = old[i];
}
void copy(double** old, double** newmatrix, const int n)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			newmatrix[i][j] = old[i][j];
}
void fulling(double** matrix, const int n)
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			std::cin >> matrix[i][j];
}
void fulling(double* vector, const int n)
{
	for (int i = 0; i < n; ++i)
		std::cin >> vector[i];
}
void print(double** matrix, const int n)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			std::cout << std::setw(12) << std::left << matrix[i][j];
		}
		std::cout << std::endl;
	}
}
void print(double* vector, const int n)
{
	for (int i = 0; i < n; ++i)
		std::cout << vector[i] << std::endl;
}

int gauss(double** A, double* B, double* X, const int n)
{
	for (int k = 0; k < n; ++k)
	{
		if (!ChangeOfStroks(A, B, n, k))
		{
			std::cout << "Matrix is degenerate" << std::endl;
			return -1;
		};
		division(A, B, n, k);
		DoZero(A, B, n, k);
	}


	for (int k = n - 1; k >= 0; --k)
	{
		X[k] = B[k];
		for (int i = 0; i < k; ++i)
			B[i] = B[i] - A[i][k] * X[k];
	}

	return 0;
}

bool ChangeOfStroks(double **A, double *B, const int n, int k)
{
	double MAX = abs(A[k][k]);
	int index = k;
	int i = k + 1;
	for (; i < n; ++i)
	{
		if (abs(A[i][k]) > MAX)
		{
			MAX = abs(A[i][k]);
			index = i;
		}
	}
	if (MAX)
	{
		if (index != k)
		{
			std::swap(A[index], A[k]);
			std::swap(B[index], B[k]);
		}
		return true;
	}
	else
		return false;
}
void division(double** A, double* B, const int n, int IndexNow)
{
	double del = A[IndexNow][IndexNow];
	for (int i = IndexNow; i < n; i++)
		A[IndexNow][i] /= del;
	B[IndexNow] /= del;
}
void DoZero(double** A, double* B, const int n, int IndexNow)
{
	for (int i = IndexNow + 1; i < n; ++i)
	{
		double M = A[i][IndexNow];
		for (int j = IndexNow; j < n; j++)
		{
			A[i][j] -= M * A[IndexNow][j];
		}
		B[i] -= M * B[IndexNow];
	}
}

void calculateError(double** A, double* X, double* B, const int n)
{
	double* rez = new double[n];
	residualVector(A, X, B, rez, n);
	std::cout << "Rezidual vector is:" << std::endl;
	print(rez, n);
	std::cout << "Norm is: " << findNorm(rez, n) << std::endl;

	AX(A, X, rez, n);
	double* XNew = new double[n];
	gauss(A, rez, XNew, n);
	std::cout << "Relative error is:" << findRelative(X, XNew, n) << std::endl;
	delete[] XNew;
}
void residualVector(double** A, double* X, double* B, double* rez, const int n)
{
	AX(A, X, rez, n);
	for (int i = 0; i < n; i++)
		rez[i] -= B[i];
}
void AX(double** A, double* X, double* rez, const int n)
{
	for (int i = 0; i < n; ++i)
	{
		rez[i] = 0;
		for (int j = 0; j < n; ++j)
			rez[i] += A[i][j] * X[j];
	}
}
double findNorm(double* arr, const int n)
{
	double MAX = 0;
	for (int i = 0; i < n; ++i)
		if (MAX < abs(arr[i]))
			MAX = abs(arr[i]);
	return MAX;
}
double findRelative(double* X, double* XNew, const int n)
{
	double MAX1 = 0;
	for (int i = 0; i < n; i++)
		if (MAX1 < abs(X[i] - XNew[i]))
			MAX1 = abs(X[i] - XNew[i]);
	double MAX2 = 0;
	for (int i = 0; i < n; i++)
		if (MAX2 < abs(X[i]))
			MAX2 = abs(X[i]);
	return MAX1 / MAX2;
}
