#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Gauss.cpp"
#include "Eiler.h"

double calculateFunction(double *u, double t, const int index)
{
	switch (index)
	{
	case 0:
		return -u[0] * u[1] + ((t < 1e-9) ? 0.0 : (sin(t) / t));
	case 1:
		return -u[1] * u[1] + (3.125 * t) / (1 + t * t);
	default:
		return 0.0;
	}
}

void eilerExplicit(double *u, const int n)
{
	double eps = 1e-3;
	double tau = 0, tauMax = 1e-2;
	double t = 0, T = 1;

	double *y = createVector(n);
	copy(u, y, n);

	double *fTemp = createVector(n);
	int k = 0;
	do
	{
		for (int i = 0; i < n; ++i)
			fTemp[i] = calculateFunction(y, t, i);

		if (eps / (std::abs(fTemp[0]) + eps / tauMax) < eps / (std::abs(fTemp[1]) + eps / tauMax))
			tau = eps / (std::abs(fTemp[0]) + eps / tauMax);
		else
			tau = eps / (std::abs(fTemp[1]) + eps / tauMax);

		for (int i = 0; i < n; ++i)
			y[i] += tau * std::abs(fTemp[i]);

		t += tau;

		for (int i = 0; i < n; i++)
		{
			std::cout << "y" << i + 1 << " = " << y[i] << std::setw(10);
		}
		std::cout << "t" << k << " = " << t << std::endl;

		k++;

	} while (t < T);

	deleteVector(fTemp);
	deleteVector(y);
}
