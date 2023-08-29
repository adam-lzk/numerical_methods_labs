#pragma once
#include<iostream>
#include<iomanip>

double** createMatrix(int n);
void deleteMatrix(double **matrix, int n);
double* createVector(int n);
void deleteVector(double* vector);
int gauss(double**, double*, double*, const int);
void fulling(double**, const int);
void print(double**, const int);
void fulling(double*, const int);
void print(double*, const int);
bool ChangeOfStroks(double**, double*, const int, int);
void division(double**, double*, const int, int);
void DoZero(double**, double*, const int, int);
void residualVector(double** A, double* X, double* B, double* rez, const int n);
double findNorm(double* arr, const int n);
void calculateError(double** A, double* X, double* B, const int n);
void copy(double* old, double* newvec, const int n);
void copy(double** old, double** newmatrix, const int n);
void AX(double** A, double* X, double* rez, const int n);
double findRelative(double* X, double* XNew, const int n);
