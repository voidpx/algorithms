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


void printmn(double *a, int m, int n) {
	printf("[");
	for (int i = 0; i < m; ++i) {
		printf("[");
		int ri = i * n;
		for (int j = 0; j < n; ++j) {
			double f = a[ri+j];
			char *s = j == n - 1 ? "]" : ", ";
			printf("%.6lf%s", f, s);
		}
		if (i != m-1) {
			printf(",\n");
		}
	}
	printf("]\n");
}

void print(double *m, int n) {
	printmn(m, n, n);
}

/*
 * m1: m x n, m2: n * k, o: m * k
 */
void matmul(double *m1, double *m2, double *o, int m, int n, int k) {
	double v;
	for (int i = 0; i < m;++i) {
		for (int j = 0; j < k; ++j) {
			v = 0.0;
			for (int x = 0; x < n; ++x) {
				v += m1[i * n + x] * m2[x * k + j];
			}
			o[i * k + j] = v;
		}
	}
}

void matmuln(double *m1, double *m2, double *o, int n) {
	matmul(m1, m2, o, n, n, n);
}

void transpose(double *a, double *o, int m, int n) {
	for (int i =0;i<m;i++) {
		for (int j = 0; j< n;j++) {
			o[j*m + i] = a[i * n + j];
		}
	}
}

void inverse(double *m, double *o, int n) {
	am a;
	a.n = n;
	memcpy(o, m, n * n * sizeof(double));
	a.m = o;
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
	memcpy(o, a.I, sizeof(double) * n * n);
	free(a.I);
	free(a.buf);
}

void cprod3d(double *m, double *n, double *o) {
	o[0] = m[1]*n[2]-n[1]*m[2];
	o[1] = m[2]*n[0]-n[2]*m[0];
	o[2] = m[0]*n[1]-n[0]*m[1];
}

double dprod(double *v1, double *v2, int n) {
	double v = 0;
	for (int i = 0; i < n; ++i) {
		v+=v1[i]*v2[i];
	}
	return v;
}

void scale(double *v, double f, int n) {
	for (int i=0;i<n;++i) {
		v[i] *= f;
	}
}

void rotate2d(double *v, double angle, double *o) {
	double a = angle * M_PI/180;
	double s = sin(a), c=cos(a);
	o[0] = c*v[0] - s*v[1];
	o[1] = s*v[0] + c*v[1];
}

void rotate3d(double *v, double *k, double angle, double *o) {
	memcpy(o, v, 3 * sizeof(double));
	double a = angle * M_PI / 180;
	double kx=k[0],ky=k[1],kz=k[2];
	double len = sqrt(kx*kx + ky*ky +kz*kz);
	kx/=len, ky/=len, kz/=len;
	double nk[3];
	nk[0] = kx, nk[1] = ky, nk[2]=kz;
	double s=sin(a), c=cos(a);
	scale(o, c, 3);
	double kxv[3];
	cprod3d(nk, v, kxv);
	scale(kxv, s, 3);
	double kvc = dprod(nk, v, 3) * (1-c);
	scale(nk, kvc, 3);
	o[0] = o[0]+kxv[0]+nk[0];
	o[1] = o[1]+kxv[1]+nk[1];
	o[2] = o[2]+kxv[2]+nk[2];
}

double det(double *m, int n) {
	if (n==1) {
		return m[0];
	}
	if (n==2) {
		return m[0]*m[3] - m[1]*m[2];
	}
	// XXX: avoid allcation in recursion
	double *sub = malloc((n-1)*(n-1)*sizeof(double));
	double d = 0;
	for (int i=0;i<n;++i) {
		double *t = sub;
		for (int j=0;j<n;++j) {
			if (j!=i) {
				memcpy(t, &m[j*n+1], (n-1)*sizeof(double));
				t+=(n-1);
			}
		}
		double s=m[i*n]*det(sub, n-1) ;
		d+= (i%2 ? -s : s);
	}
	free(sub);
	return d;
}

/*
 * least square fit. a: m x n, b: m x 1, o: n x 1.
 *
 */
void lsf(double *a, double *b, int m, double *o, int n) {
	assert(m >= n);
	double *at = malloc(m * n * sizeof(double));
	assert(at);
	transpose(a, at, m, n);
	double *ata = malloc(n * n * sizeof(double));
	assert(ata);
	matmul(at, a, ata, n, m, n);
	double *atai = malloc(n * n * sizeof(double));
	assert(atai);
	inverse(ata, atai, n);
	double *atb = malloc(n * sizeof(double));
	matmul(at, b, atb, n, m, 1);
	matmul(atai, atb, o, n, n, 1);
	free(atb);
	free(atai);
	free(ata);
	free(at);
}

int main() {
	double a[16] = {1,2,3,4,5,6,7,8,0,99,33,2,11,23,76,9};
	print(a, 3);
	double d = det(a,3);
	printf("%lf\n", d);

	double v[3]={1,1,1};
	double k[3]={0,0,1};
	double o[3];
	rotate3d(v, k, 45, o);

	printmn(o, 3, 1);

	double v2[2]={1,1};
	double r2[2];
	rotate2d(v2, 45, r2);

	printmn(r2, 2, 1);




#define M 4
#define N 3
	double m[M * N] = {1., 2., 3., 4., 5., 6., 12, 9, 3, 44, 89, 11};
//	double inv[N * N];
	printmn(m, M, N);
//	inverse(m, inv, N);
//	print(inv, N);
	double t[N];
//	transpose(m, t, M, N);
//	matmuln(m, inv, t, N);

	double b[M] = {1,2,3,4};
//	printmn(b, M, 1);

	double y[N];
	lsf(m, b, M, y, N);

	printmn(y, N, 1);
	return 0;
}
