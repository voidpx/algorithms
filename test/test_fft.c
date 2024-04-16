#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "../src/fft.h"
#include "minunit.h"

static void print(char* head, double *v, uint32_t len) {
	printf("%s:\n", head);
	for (uint32_t i = 0; i < len; ++i) {
		printf("%f%s", v[i], i==len-1 ? "" : ", ");
	}
	printf("\n");
}

static void print_comp(char*head, complex *v, uint32_t len) {
	printf("%s:\n", head);
	for (uint32_t i = 0; i < len; ++i) {
		printf("%f%+fj%s", creal(v[i]), cimag(v[i]), i==len-1 ? "" : ", ");
	}
	printf("\n");
}

static void test_fft() {
#define N 1000
	double data[N];
	double s = 1.0/N;

	for (int i = 0; i < N; ++i) {
		double t = i*s;
		double y = 3 * sin(2 * M_PI * t);
		y += sin(2 * M_PI * 4 * t);
		y += 0.5 * sin(2 * M_PI * 7 * t);
		data[i] = y;
	}

//	print("data", data, N);

	FFT p = dft(data, N);
//	print_comp("DFT", p.value, p.len);


	FFT pf = fft(data, N);
//	print_comp("FFT", pf.value, pf.len);

	assert(p.len == pf.len);
	for (int i = 0; i < p.len; ++i) {
		assert(abs(creal(p.value[i]) - creal(pf.value[i])) < 10e-8 && abs(cimag(p.value[i]) - cimag(pf.value[i])) < 10e-8);
	}

	double *d = ifft(pf);

//	print("IFFT:", d, N);

	for (int i = 0; i < N; ++i) {
		assert(abs(d[i] - data[i]) < 10e-8);
	}

	free(d);
	free(p.value);
	free(pf.value);
}

int main() {
	test_fft();
}
