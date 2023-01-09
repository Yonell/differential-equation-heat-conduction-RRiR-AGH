#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
//#define ARMA_USE_SUPERLU
#include <armadillo>
#include <conio.h>

double xk_func(double max, double min, int k, int n) {
	return (max - min) * k / n + min;
}

double ek(double x, double max, double min, int k, int n) {
	double xk_minus_1 = xk_func(max, min, k - 1, n);
	double xk = xk_func(max, min, k, n);
	double xk_plus_1 = xk_func(max, min, k + 1, n);
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

double ek_derivative(double x, double max, double min, int k, int n) {
	double xk_minus_1 = xk_func(max, min, k - 1, n);
	double xk = xk_func(max, min, k, n);
	double xk_plus_1 = xk_func(max, min, k + 1, n);
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

double Lj(int j, int n) {
	return j==0 ? -20 * ek(0, 2, 0, j, n) : 0;
}

double trapeze_method_for_the_product_two_ek_derivatives(int i, int j, double a, double b, int n, double step) {
	double sum = 0;
	for (double x = a; x <= b; x += step) {
		sum += (ek_derivative(x, 2, 0, i, n) * ek_derivative(x, 2, 0, j, n) 
			+ ek_derivative(x + step, 2, 0, i, n) * ek_derivative(x + step, 2, 0, j, n)) 
			* step / 2;
	}
	return sum;
}

double Bij(int i, int j, int n) {
	return (abs(i - j) < 2) ? (
		- ek(1, 2, 0, j, n) * ek_derivative(1, 2, 0, i, n)
		- ek(0, 2, 0, j, n) * ek(0, 2, 0, i, n)
		- 2 * ek(2, 2, 0, j, n) * ek_derivative(2, 2, 0, i, n)
		+ 2 * ek(1, 2, 0, j, n) * ek_derivative(1, 2, 0, i, n)
		+ trapeze_method_for_the_product_two_ek_derivatives(i, j, 0, 1, n, 1 / (10.0 * n))
		+ 2 * trapeze_method_for_the_product_two_ek_derivatives(i, j, 1, 2, n, 1 / (10.0 * n))
	) : 0;
}

double u(double x, arma::vec w, int n) {
	double val = 0;
	arma::vec e_values(n);
	for (int i = 0; i < n; i++) {
		e_values(i) = w(i)*ek(x, 2, 0, i, n);
	}
	return arma::sum(e_values);
}

void show_help() {
	std::cout << "\nThis is help for the heat conduction program made by Filip Piskorski.\n";
	std::cout << "You can run this program with or withour argument.\n\n";
	std::cout << "Available options: \n\n";
	std::cout << "  -m / --show_matrices            - Display the matricies.\n";
	std::cout << "  -a / --show_answer              - Display the W matrix.\n";
	std::cout << "  -g / --show_graph               - Show graph (requires gnuplot).\n";
	std::cout << "  -sp / --use_sparse_algorithm    - Use sparse algorith while solving (Doesn't work :/ ).\n";
	std::cout << "  -n / --number_of_elements       - Specify the number of elements before starting the program.\n";
	std::cout << "  -? / -h / --help                - Show help.\n";

	std::cout << "\n";

	return;
}

int main(int argc, char* argv[]) {
	int n_elements;
	bool show_matrices = false;
	bool show_answer = false;
	bool show_graph = false;
	bool use_sparse = false;
	bool use_argument_n = false;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--show_matrices") == 0)
			show_matrices = true;
		if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--show_answer") == 0)
			show_answer = true;
		if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--show_graph") == 0)
			show_graph = true;
		if (strcmp(argv[i], "-sp") == 0 || strcmp(argv[i], "--use_sparse_algorithm") == 0){
			use_sparse = true;
			std::cout << "Sparse matricies are not available at the moment, because for some reason superlu library doesn't work.\n";
			return -1;
		}
		if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--number_of_elements") == 0) {
			use_argument_n = true;
			n_elements = strtol(argv[i + 1], nullptr, 10);
			i++;
		}
		if (strcmp(argv[i], "-?") == 0 || strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			show_help();
			return 0;
		}

	}
	if (!use_argument_n) {
		std::cout << "Enter the number of elements: ";
		std::cin >> n_elements;
		std::cout << "\n";
	}
	arma::mat B(n_elements, n_elements);
	arma::vec L(n_elements);
	arma::vec W;
	for (int j = 0; j < n_elements; j++)
	{
		for (int i = 0; i < n_elements; i++) {
			B(i,j) = Bij(i, j, n_elements);
		}
		L(j) = Lj(j, n_elements);
	}

	if (use_sparse) {
		//W = arma::spsolve(arma::sp_mat(B), L);
	}
	else {
		W = arma::solve(B, L);
	}

	std::cout << "\n Solving equation: BW = L \n\n";

	if (show_matrices) {
		std::cout << "B =\n" << B << "\n";
		std::cout << "L =\n" << L << "\n";
	}
	if (show_answer) {
		std::cout << "W =\n" << W;
	}

	if (show_graph) {
		FILE* temp = NULL;
		FILE* gnupipe = _popen("gnuplot -persistent", "w");
		temp = fopen("data.tmp", "w");
		fprintf(temp, "%f %f\n", 0, u(0, W, n_elements));
		for (double i = 0.01; i <= 2; i += 0.01) {
			fprintf(temp, "%f %f\n", i, u(i, W, n_elements));
		}
		fprintf(temp, "%f %f\n", 2.0, 0.0);
		fprintf(gnupipe, "%s\n", "plot 'data.tmp' with lines lc rgb \"#ff0000\"");
	}
	return 0;
}