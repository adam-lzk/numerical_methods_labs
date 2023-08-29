#pragma once
#include<cmath>
#include"Gauss.h"

typedef double(*pf)(double*, double*, double, double);

bool Newton(double *yk1, double *yk, double tk, double Tau, const int n);
double f1(double *uk1, double *uk, double t, double Tau);
double f2(double *uk1, double *uk, double t, double Tau);
double Differential(pf f, double *uk1, double *uk, double t, double Tau, int n);
void inversionVector(double* vector, const int n);
void updateX(double* X, double* dX, const int n);
double changeOfUpdateDX(double* X, double* dX, const int n);
