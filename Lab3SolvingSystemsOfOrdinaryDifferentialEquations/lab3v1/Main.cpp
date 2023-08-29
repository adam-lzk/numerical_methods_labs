#pragma once
#include"Gauss.cpp"
#include"eilerExplicit.cpp"
#include"eilerImplicit.cpp"

int main()
{
	int n = 2;
	double* u = createVector(n);
	u[0] = 0;
	u[1] = -0.412;

	eilerImplicit(u, n);
	std::cout << "==========================" << std::endl << "==========================" <<std::endl;
	eilerExplicit(u, n);
	return 0;
}
