#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
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

static double *test_data(uint32_t n) {
	double *data = malloc(sizeof(double) * n);
	double s = 1.0/n;

	for (int i = 0; i < n; ++i) {
		double t = i*s;
		double y = 3 * sin(2 * M_PI * t);
		y += sin(2 * M_PI * 4 * t);
		y += 0.5 * sin(2 * M_PI * 7 * t);
		data[i] = y;
	}
	return data;
}

static void test_fft(double *data, uint32_t n, FFT p) {
//	print("data", data, n);
//	print_comp("DFT", p.value, p.len);

	FFT pf = fft(data, n);
//	print_comp("FFT", pf.value, pf.len);

	assert(p.len == pf.len);
	for (int i = 0; i < p.len; ++i) {
		assert(abs(creal(p.value[i]) - creal(pf.value[i])) < 10e-8 && abs(cimag(p.value[i]) - cimag(pf.value[i])) < 10e-8);
	}

	double *d = ifft(pf);

//	print("IFFT:", d, n);

	for (int i = 0; i < n; ++i) {
		assert(abs(d[i] - data[i]) < 10e-8);
	}

	free(d);
	free(pf.value);
}

typedef union TR {
	FFT fft;
	double * real;
} TR;

typedef struct func {
	TR (*f)(va_list);
} func;

static TR time_it(char *s, func f, ...) {
	struct timespec t0, t1;
	clock_gettime(CLOCK_MONOTONIC_RAW, &t0);
	va_list va;
	va_start(va, f);
	TR r = f.f(va);
	va_end(va);
	clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
	double t = (t1.tv_sec - t0.tv_sec) + ((double)(t1.tv_nsec - t0.tv_nsec))/1.0e+9;
	printf("%s, time: %.6fs\n", s, t);
	return r;
}

static FFT tf(FFT (*f)(double *, uint32_t), va_list a) {
	double *data = va_arg(a, double*);
	uint32_t n =  va_arg(a, uint32_t);
	return f(data, n);
}

static TR dftf(va_list a) {
	return (TR){tf(dft, a)};
}

static TR fftf(va_list a) {
	return (TR){tf(fft, a)};
}

int main() {
#define N 2000
	struct timespec t0, t1;
	double *data = test_data(N);
	func df = {dftf};
	TR dftr = time_it("DFT", df, data, (uint32_t)N);

	// check correctness
	test_fft(data, N, dftr.fft);

	func ff = {fftf};
	TR fftr = time_it("FFT", ff, data, (uint32_t)N);

	free(dftr.fft.value);
	free(fftr.fft.value);

}
