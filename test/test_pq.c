#include <stdio.h>
#include <string.h>
#include "../src/priorityqueue.h"
#include "minunit.h"

int comp(void *a, void *b) {
	return *((int*) a) - *((int*) b);
}

#define ARR_LEN(a) sizeof(a)/sizeof(a[0])

static void test_pq() {
	PriorityQueue *p = new_pq(10, comp);
	int a[] = { 2, 3, 11, 0, 55, 90, 23, 45, 4, 100, 3, 4, 6, 7 };

	for (int i = 0; i < ARR_LEN(a); i++) {
		enqueue(p, &a[i]);
	}
	void *o;
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 0);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 2);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 3);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 3);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 4);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 4);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 6);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 7);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 11);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 23);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 45);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 55);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 90);
	o = dequeue(p);
	mu_assert("failed", *((int* )o) == 100);
	o = dequeue(p);
		mu_assert("failed", o == NULL);
	free_pq(p);
}

int main(void) {
	test_pq();
}
