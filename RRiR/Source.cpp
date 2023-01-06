#include <iostream>
#include <string>

float xk_func(float max, float min, int k, int n) {
	return (max - min) * k / n + min;
}

float ek(float x, float max, float min, int k, int n) {
	float xk_minus_1 = xk_func(max, min, k - 1, n);
	float xk = xk_func(max, min, k, n);
	float xk_plus_1 = xk_func(max, min, k = 1, n);
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
	float xk_plus_1 = xk_func(max, min, k = 1, n);
	if (x > xk_minus_1 && x < xk) {
		return 1 / (xk - xk_minus_1);
	}
	else if (x > xk && x < xk_plus_1) {
		return -1 / (xk_plus_1 - xk);
	}
	else {
		return 0;
	}
}

int main(std::string args[]) {
	
}