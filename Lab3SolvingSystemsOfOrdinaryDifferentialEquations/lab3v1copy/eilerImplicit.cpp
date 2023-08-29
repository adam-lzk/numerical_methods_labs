#pragma once
#include <iostream>
#include <iomanip>
#include "Newton.cpp"
#include "Eiler.h"

void eilerImplicit(double* u, const int n)
{
	double Tauk, Tauk_1, Tauk1;
	double eps = 1e-3, T = 1, TauMax = 1e-2, TauMin = 1e-2;
	double tk = 0, tk1 = 0;
	double* epsk = createVector(n);
	double* yk = createVector(n);
	double* yk_1 = createVector(n);
	double* yk1 = createVector(n);

	for (int i = 0; i < n; i++)
		yk[i] = yk_1[i] = yk1[i] = u[i];

	Tauk = Tauk_1 = TauMin;

	int k = 1;
	while (tk < T)
	{
		tk1 = tk + Tauk;
		Newton(yk1, yk, tk, Tauk, n);

		for (int i = 0; i < n; ++i)
			epsk[i] = -(Tauk / (Tauk + Tauk_1))*(yk1[i] - yk[i] - Tauk * (yk[i] - yk_1[i]) / Tauk_1);

		bool flag = false;
		for (int i = 0; i < n && !flag; i++)
			if (std::fabs(epsk[i]) > eps && !flag)
			{
				Tauk /= 2;
				tk1 = tk;
				for (int j = 0; j < n; j++)
					yk1[j] = yk[j];
				flag = true;
			}
		if (flag)
			continue;
		double *temp = createVector(n);

		for (int i = 0; i < n; i++)
		{
			temp[i] = sqrt(eps / std::abs(epsk[i]))*Tauk;
		}
		for (int i = 0; i < n - 1; i++)
		{
			if (temp[i] > temp[i + 1])
				Tauk1 = temp[i + 1];
			else
				Tauk1 = temp[i];
		}
		deleteVector(temp);

		if (Tauk1 > TauMax)
			Tauk1 = TauMax;



		std::cout << "y" << k + 1 << "[0] = " << yk1[0]
			<< std::setw(10) << "y" << k + 1 << "[1] = " << yk1[1]
			<< std::setw(10) << "t" << k + 1 << " = " << tk1 << std::endl;


		for (int i = 0; i < n; i++)
		{
			yk_1[i] = yk[i];
			yk[i] = yk1[i];
		}
		Tauk_1 = Tauk;
		Tauk = Tauk1;
		tk = tk1;

		k++;
	}

	deleteVector(epsk);
	deleteVector(yk);
	deleteVector(yk1);
	deleteVector(yk_1);
}
