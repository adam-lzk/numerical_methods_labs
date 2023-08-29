#pragma once
#include"Newton.h"
#include"Gauss.cpp"

bool Newton(double *yk1, double *yk, double tk, double Tau, const int n)
{
	double d1, d2;
	const double eps = 1e-9;
	d1 = d2 = eps;
	const int NIT = 100;

	double** Jacobi = createMatrix(n * n);
	double* rezidual = createVector(n);
	double* dyk1 = createVector(n);

	int k = 0;
	for (; k < NIT; k++)
	{
		rezidual[0] = f1(yk1, yk, tk, Tau);
		rezidual[1] = f2(yk1, yk, tk, Tau);

		Jacobi[0][0] = Differential(f1, yk1, yk, tk, Tau, 1);
		Jacobi[0][1] = Differential(f1, yk1, yk, tk, Tau, 2);
		Jacobi[1][0] = Differential(f2, yk1, yk, tk, Tau, 1);
		Jacobi[1][1] = Differential(f2, yk1, yk, tk, Tau, 2);

		inversionVector(rezidual, n);
		if (gauss(Jacobi, rezidual, dyk1, n) == -1)
			return false;

		if (findNorm(rezidual, n) < d1 && changeOfUpdateDX(yk1, dyk1, n) < d2)
			break;
		updateX(yk1, dyk1, n);

	}
	if (k == 100)
		std::cout << "IER = 2\n";


	deleteMatrix(Jacobi, n);
	deleteVector(rezidual);
	deleteVector(dyk1);
	return true;

}


double f1(double *uk1, double *uk, double t, double Tau)
{
	return uk1[0] - uk[0] - Tau * (-uk1[0] * uk1[1] + ((t < 1e-9) ? 0.0 : (sin(t) / t)));
}
double f2(double *uk1, double *uk, double t, double Tau)
{
	return uk1[1] - uk[1] - Tau * (-uk1[1] * uk1[1] + (3.125*t) / (1 + t * t));
}

typedef double(*pf)(double*, double*, double, double);

double Differential(pf f, double *uk1, double *uk, double t, double Tau, int n)
{
	double dx = 1e-9;
	double* D = new double[3];
	for (int i = 0; i < 2; i++)
	{
		if (i == n - 1)
		{
			D[i] = uk1[i] + dx;
			i++;
		}
		D[i] = uk1[i];
	}

	double F = f(uk1, uk, t, Tau);
	double dF = f(D, uk, t, Tau);

	delete[] D;
	return (dF - F) / dx;
}

void inversionVector(double* vector, const int n)
{
	for (int i = 0; i < n; i++)
		vector[i] *= -1;
}

void updateX(double* X, double* dX, const int n)
{
	for (int i = 0; i < n; i++)
		X[i] += dX[i];
}

double changeOfUpdateDX(double* X, double* dX, const int n)
{
	double norm;
	double max = 0;
	double* temp = createVector(n);
	double* rezidualOfXk = createVector(n);
	double* vec = createVector(n);
	for (int i = 0; i < n; i++)
		temp[i] = X[i] + dX[i];
	for (int i = 0; i < n; i++)
		if (abs(temp[i] < 1))
			vec[i] = dX[i];
		else
		{
			rezidualOfXk[i] = dX[i];
			vec[i] = rezidualOfXk[i] /= temp[i];
		}
	norm = findNorm(vec, n);
	deleteVector(vec);
	deleteVector(rezidualOfXk);
	deleteVector(temp);
	return norm;
}
