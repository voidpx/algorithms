#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fft.h"

// a simple/naive implementation of the FFT algorithm:
// https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
// normalization factor is 1/N, applied to the forward pass

//https://graphics.stanford.edu/%7Eseander/bithacks.html#RoundUpPowerOf2
static uint32_t po2(uint32_t n) {
	n--;
	n |= n>>1;
	n |= n>>2;
	n |= n>>4;
	n |= n>>8;
	n |= n>>16;
	n++;
	return n;
}

typedef struct V {
	complex *v;
	uint32_t len;
} V;

static complex _v(V* v, uint32_t index) {
	if (index >= v->len) {
		return 0.0;
	}
	return v->v[index];
}

typedef struct C {
	complex *p;
	uint32_t n;
} C;

/*
 * dir, -1: forward FFT, 1: inverse FFT
 */
static C _c(uint32_t po2len, int dir) {
	// XXX: could be further improved to only calculate pi/4
	int n = po2len / 4;
	double a = dir * M_PI_2 / n;
	complex *d = malloc(sizeof(complex) * po2len);
	for (int i = 0; i < po2len; ++i) {
		int q = i / n;
		int r = i % n;
		switch (q) {
		case 0:
			double s = i * a;
			d[i] = cos(s) + sin(s)*I;
			continue;
		case 1:
			d[i] = -dir*cimag(d[r]) + dir*creal(d[r]) * I;
			break;
		case 2:
			d[i] = -creal(d[r]) - cimag(d[r]) * I;
			break;
		case 3:
			d[i] = dir*cimag(d[r]) - dir*creal(d[r]) * I;
			break;
		}
	}
	return (C){d, po2len};
}

// to confirm the above was correct
//static C _c_test(uint32_t po2len) {
//	complex *d = malloc(sizeof(complex) * po2len);
//	double a =2*M_PI / po2len;
//	for (int i = 0; i < po2len; ++i) {
//		d[i] = cexp(-I*2*M_PI * i / po2len);
//	}
//	return (C ) { d, po2len };
//}

/*
 * recursively calculate the fft. NOTE: subtle difference: not exactly fft in each step,
 * the odd part is not multiplied by the factor exp(-2pi*i*k/N), since in the pure recursive
 * version, it's just the odd part is rotated clockwise first and then counter-clockwise by the
 * same factor, which is not needed here.
 * @k - frequency
 * @s - start index of this round
 * @m - multiplier, power of 2
 */

typedef struct R {
	complex rk;
	complex rkn2;
} R;

/*
 * calculate the DFT for frequency k.
 */
static R _fftr(V *v, C *c, uint32_t k, uint32_t s, uint32_t m) {
	if (m >= c->n) {
		return (R){_v(v, s), 0};
	}
	uint32_t m2 = m << 1;
	R ek = _fftr(v, c, k, s, m2);
	R ok = _fftr(v, c, k, s + m, m2);
	complex t = ok.rk * c->p[(m * k) % c->n]; // twiddle factor
	return (R){ek.rk + t, m == 1 ? ek.rk - t : 0}; // only need first iteration for k + N/2
}

static complex *_fftn(V *v, uint32_t po2len) {
	if (!po2len) {
		return NULL;
	}
	assert((po2len & po2len - 1) == 0);
	C c = _c(po2len, -1);
	if (!c.p) {
		return NULL;
	}
	complex *p = calloc(po2len, sizeof(complex));
	if (!p) {
		return NULL;
	}
	for (int i=0; i < po2len/2; ++i) {
		complex even = {0};
		complex odd = {0};
		for (int e = 0, o = 1; o < c.n; e+=2,o+=2) {
			even += e==0? v->v[0] : (i == 0 ? v->v[e] : v->v[e] * c.p[(e * i) % po2len]);
			odd += i==0? v->v[o] : v->v[o] * c.p[(o * i)  % po2len];
		}
		p[i] = (even + odd)/c.n;
		p[i+po2len/2] = (even - odd)/c.n;
	}
	free(c.p);
	return p;
}

/*
 * plain DFT.
 */
static complex *_dft(V *v, uint32_t po2len) {
	if (!po2len) {
		return NULL;
	}
	assert((po2len & po2len - 1) == 0);
	C c = _c(po2len, -1);
	if (!c.p) {
		return NULL;
	}
	complex *p = malloc(sizeof(complex) * po2len);
	if (!p) {
		return NULL;
	}
	for (int i = 0; i < po2len; ++i) {
		complex s = _v(v, 0);
		for (int j = 1; j < po2len; ++j) {
			s += _v(v, j) * c.p[(j * i) % po2len];
		}
		p[i] = s /= po2len;
	}
	free(c.p);
	return p;
}

static complex* _fft(V *v, uint32_t po2len, int dir) {
	if (!po2len) {
		return NULL;
	}
	assert((po2len & po2len - 1) == 0);
	C c = _c(po2len, dir);
	if (!c.p) {
		return NULL;
	}
	complex *p = calloc(po2len, sizeof(complex));
	if (!p) {
		return NULL;
	}
	for (int i=0; i < po2len/2; ++i) {
		R ak = _fftr(v, &c, i, 0, 1);
		p[i] = dir == -1 ? ak.rk / po2len : ak.rk;
		p[i + po2len/2] = dir == -1 ? ak.rkn2 / po2len : ak.rkn2;
	}
	free(c.p);
	return p;
}

static complex *_complexify(double *v, uint32_t len) {
	complex *a = malloc(sizeof(complex) * len);
	if (!a) {
		return NULL;
	}
	for (int i = 0; i < len; ++i) {
		a[i] = v[i] + 0*I;
	}
	return a;
}

/*
 * direct implementation of DFT
 *
 */
FFT dft(double *v, uint32_t len) {
	complex *a = _complexify(v, len);
	if (!a) {
		return (FFT){0};
	}
	V val = {a, len};
	uint32_t po2len = po2(len);
	complex *r = _dft(&val, po2len);
	if (!r) {
		return (FFT){0};
	}
	free(a);
	return (FFT){r, po2len};
}

/*
 * FFT
 */
FFT fftn(double *v, uint32_t len) {
	complex *a = _complexify(v, len);
	if (!a) {
		return (FFT){0};
	}
	V val = {a, len};
	uint32_t po2len = po2(len);
	complex *r = _fftn(&val, po2len);
	if (!r) {
		return (FFT){0};
	}
	free(a);
	return (FFT){r, po2len};
}

/*
 * FFT
 */
FFT fft(double *v, uint32_t len) {
	complex *a = _complexify(v, len);
	if (!a) {
		return (FFT){0};
	}
	V val = {a, len};
	uint32_t po2len = po2(len);
	complex *r = _fft(&val, po2len, -1);
	if (!r) {
		return (FFT){0};
	}
	free(a);
	return (FFT){r, po2len};
}

/*
 * inverse FFT
 */
double *ifft(FFT c) {
	V v={c.value, c.len};
	int po2len = po2(c.len);
	complex *inv = _fft(&v, po2len, 1);
	if (!inv) {
		return NULL;
	}
	double *d = malloc(sizeof(double) * po2len);
	if (!d) {
		return NULL;
	}
	for (int i = 0; i < po2len; ++i) {
		d[i] = creal(inv[i]);
	}
	free(inv);
	return d;
}



