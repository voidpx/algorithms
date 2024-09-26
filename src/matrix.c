#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define PREC 10E-10

typedef struct aug_matrix {
	double *m;
	double *I;
	int n;
	double *buf;
} am;

void op_mul(am *m, int r, double f) {
	double *mp = &m->m[r * m->n];
	double *ip = &m->I[r * m->n];
	for (int i = 0; i < m->n; ++i) {
		mp[i] *= f;
		ip[i] *= f;
	}
}

void op_mul_add(am *m, int r1, int r2, double f) {
	double *r1p = &m->m[r1 * m->n];
	double *r2p = &m->m[r2 * m->n];
	double *r1ip = &m->I[r1 * m->n];
	double *r2ip = &m->I[r2 * m->n];
	for (int i = 0; i < m->n; ++i) {
		r2p[i] += r1p[i] * f;
		r2ip[i] += r1ip[i] * f;
	}
}

void op_swap(am *m, int r1, int r2) {
	double *p = m->buf;
	assert(p);
	size_t rl = sizeof(double) * m->n;
	int r1i = r1 * m->n, r2i = r2 * m->n;
	memcpy(p, &m->m[r1i], rl);
	memcpy(&m->m[r1i], &m->m[r2i], rl);
	memcpy(&m->m[r2i], p, rl);

	memcpy(p, &m->I[r1i], rl);
	memcpy(&m->I[r1i], &m->I[r2i], rl);
	memcpy(&m->I[r2i], p, rl);

}

static int is_zero(double d) {
	return fabs(d) < PREC;
}


void print(double *m, int n) {
	printf("[");
	for (int i = 0; i < n ; ++i) {
		printf("[");
		int ri = i * n;
		for (int j = 0; j < n; ++j) {
			double f = m[ri+j];
			char *s = j == n - 1 ? "]" : ", ";
			printf("%.6lf%s", f, s);
		}
		if (i != n-1) {
			printf(",\n");
		}
	}
	printf("]\n");
}

void matmul(double *m1, double *m2, double *t, int n) {
	double v;
	for (int i = 0; i < n;++i) {
		for (int j = 0; j < n; ++j) {
			v = 0.0;
			for (int x = 0; x < n; ++x) {
				v += m1[i * n + x] * m2[x * n + j];
			}
			t[i * n + j] = v;
		}
	}
}

void inverse(double *m, double *t, int n) {
	am a;
	a.n = n;
	memcpy(t, m, n * n * sizeof(double));
	a.m = t;
	a.I = malloc(n * n * sizeof(double));
	a.buf = malloc(n * sizeof(double));
	assert(a.I);
	assert(a.buf);
	for (int i = 0; i < n; ++i) {
		int r = i * n;
		for (int j = 0; j < n; ++j) {
			if (i == j) {
				a.I[r + j] = 1.d;
			} else {
				a.I[r + j] = 0.0;
			}
		}
	}
	for (int z = 0; z < n; ++z) {
		if (is_zero(a.m[z * n + z])){
			for (int i = z+1; i < n; ++i) {
				if (!is_zero(a.m[i * n + z])) {
					op_swap(&a, z, i);
					break;
				}
			}
		}
		assert(!is_zero(a.m[z * n + z]));
		for (int i = z+1; i < n; ++i) {
			if (!is_zero(a.m[i * n + z])) {
				double x = -a.m[i * n + z]/a.m[z * n + z];
				op_mul_add(&a, z, i, x);
			}
		}
	}

	for (int z = n-1; z >= 0; --z) {
		for (int i = z-1; i >= 0; --i) {
			if (!is_zero(a.m[i * n + z])) {
				double x = -a.m[i * n + z]/a.m[z * n + z];
				op_mul_add(&a, z, i, x);
			}
		}
		if (fabs(a.m[z * n + z] - 1.0) > PREC) {
			op_mul(&a, z, 1.0/a.m[z * n + z]);
		}
	}
	memcpy(t, a.I, sizeof(double) * n * n);
	free(a.I);
	free(a.buf);
}

int main() {
#define N 3
	double m[N * N] = {1., 2., 3., 4., 5., 6., 1., 8., 9.};
	double inv[N * N];
	print(m, N);
	inverse(m, inv, N);
	print(inv, N);
	double t[N * N];
	matmul(m, inv, t, N);
	print(t, N);
	return 0;
}
