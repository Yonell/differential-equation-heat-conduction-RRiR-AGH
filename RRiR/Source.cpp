#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#define ARMA_USE_SUPERLU true
//#include "slu_cdefs.h"
//#include "slu_Cnames.h"
//#include "slu_dcomplex.h"
//#include "slu_ddefs.h"
//#include "slu_scomplex.h"
//#include "slu_sdefs.h"
//#include "slu_util.h"
//#include "slu_zdefs.h"
//#include "superlu_enum_consts.h"
//#include "supermatrix.h"
#include <armadillo>
#include <conio.h>

float xk_func(float max, float min, int k, int n) {
	return (max - min) * k / n + min;
}

float ek(float x, float max, float min, int k, int n) {
	float xk_minus_1 = xk_func(max, min, k - 1, n);
	float xk = xk_func(max, min, k, n);
	float xk_plus_1 = xk_func(max, min, k + 1, n);
	if (x > xk_minus_1 && x <= xk) {
		return (x - xk_minus_1) / (xk - xk_minus_1);
	}
	else if (x > xk && x < xk_plus_1) {
		return (xk_plus_1 - x) / (xk_plus_1 - xk);
	}
	else {
		return 0;
	}
}

float ek_derivative(float x, float max, float min, int k, int n) {
	float xk_minus_1 = xk_func(max, min, k - 1, n);
	float xk = xk_func(max, min, k, n);
	float xk_plus_1 = xk_func(max, min, k + 1, n);
	if (x > xk_minus_1 && x < xk) {
		return 1 / (xk - xk_minus_1);
	}
	else if (x > xk && x < xk_plus_1) {
		return -1 / (xk_plus_1 - xk);
	}
	else {
		return 0;
	}
	std::cout << "heh";
}

float Lj(int j, int n) {
	return j==0 ? -20 * ek(0, 2, 0, j, n) : 0;
}

float trapeze_method_for_the_product_two_ek_derivatives(int i, int j, float a, float b, int n, float step) {
	float sum = 0;
	for (float x = a; x <= b; x += step) {
		sum += (ek_derivative(x, 2, 0, i, n) * ek_derivative(x, 2, 0, j, n) 
			+ ek_derivative(x + step, 2, 0, i, n) * ek_derivative(x + step, 2, 0, j, n)) 
			* step / 2;
	}
	return sum;
}

float Bij(int i, int j, int n) {
	return (abs(i - j) < 2) ? (
		- ek(1, 2, 0, j, n) * ek_derivative(1, 2, 0, i, n)
		- ek(0, 2, 0, j, n) * ek(0, 2, 0, i, n)
		- 2 * ek(2, 2, 0, j, n) * ek_derivative(2, 2, 0, i, n)
		+ 2 * ek(1, 2, 0, j, n) * ek_derivative(1, 2, 0, i, n)
		+ trapeze_method_for_the_product_two_ek_derivatives(i, j, 0, 1, n, 1 / (10.0 * n))
		+ 2 * trapeze_method_for_the_product_two_ek_derivatives(i, j, 1, 2, n, 1 / (10.0 * n))
	) : 0;
}

float u(float x, arma::vec w, int n) {
	float val = 0;
	arma::vec e_values(n);
	for (int i = 0; i < n; i++) {
		e_values(i) = w(i)*ek(x, 2, 0, i, n);
	}
	return arma::sum(e_values);
}

int main(int argc, char* argv[]) {
	int n_elements;
	//if (argc == 0) {
		std::cout << "Enter the number of elements: ";
		std::cin >> n_elements;
		std::cout << "\n";
	//}
	//else if (argc == 1) {
	//	n_elements = strtol(argv[0], nullptr, 10);
	//}
	//else {
	//	std::cout << "bad argument/s!!! :c \n";
	//	return -1;
	//}
	arma::mat B(n_elements, n_elements);
	arma::vec L(n_elements);
	for (int j = 0; j < n_elements; j++)
	{
		for (int i = 0; i < n_elements; i++) {
			B(i,j) = Bij(i, j, n_elements);
		}
		L(j) = Lj(j, n_elements);
	}
	std::cout << B << "\n";
	std::cout << L << "\n";

	arma::vec Ans;
	Ans = arma::solve(B, L);
	std::cout << Ans;


	FILE* temp = NULL;
	FILE* gnupipe = _popen("gnuplot -persistent", "w");
	temp = fopen("data.tmp", "w");
	for (float i = 0; i <= 2; i += 0.01) {
		fprintf(temp, "%f %f\n", i, u(i, Ans, n_elements));
	}
	fprintf(temp, "%f %f\n", 2.0, 0.0);
	fprintf(gnupipe, "%s\n", "plot 'data.tmp'");

	std::cout << "Press any key to show the graph... ";
	_getch();
}