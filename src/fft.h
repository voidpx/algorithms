#ifndef SRC_FFT_H_
#define SRC_FFT_H_

#include <stdint.h>
#include <complex.h>

typedef struct FFT {
	complex *value;
	uint32_t len;
} FFT;

/*
 * DFT, it's up to the caller to free the pointer in the result, whose length is len rounded to the next power of 2.
 * NULL if no memory could be allocated.
 */
FFT dft(double *v, uint32_t len);

/*
 * FFT.
 */
FFT fft(double *v, uint32_t len);

/*
 * inverse FFT.
 */
double *ifft(FFT c);

#endif /* SRC_FFT_H_ */
