#include "priorityqueue.h"

#include <stdlib.h>
#include <string.h>

#define DEF_SIZE 0xF

static void _shift_up(PriorityQueue *p, size_t i, void *o);
static void _shift_down(PriorityQueue *p, size_t i, void *o);
static int _index_of(PriorityQueue *p, void *o);
static int _remove_at(PriorityQueue *p, size_t i);
static int _shrink(PriorityQueue *p);

PriorityQueue* new_pq(size_t capacity, int (*comparator)(void *, void *));
void free_pq(PriorityQueue *p);
int enqueue(PriorityQueue *p, void *obj);
int remove_element(PriorityQueue *p, void *obj);
void* dequeue(PriorityQueue *p);

PriorityQueue* new_pq(size_t capacity, int (*comparator)(void *, void *)) {
	PriorityQueue* p = malloc(sizeof (PriorityQueue));
	void* queue = malloc(capacity * sizeof (void *));
	if (!queue) {
		return NULL;
	}
	p->queue = queue;
	p->capacity = capacity;
	p->comparator = comparator;
	p->size = 0;
	return p;
}

void free_pq(PriorityQueue *p) {
	free(p->queue);
	free(p);
}

int grow(PriorityQueue *p) {
	int s;
	if (p->capacity == 0) {
		s = DEF_SIZE;
	} else {
		s = p->capacity + ((p->capacity + 1) >> 1);
	}
	p->queue = realloc(p->queue, s * sizeof(void *));
	p->capacity = s;
	return 0;
}

int enqueue(PriorityQueue *p, void* obj) {
	if (p->size >= p->capacity) {
		if (grow(p)) {
			return ERR;
		}
	}
	p->queue[p->size] = obj;
	_shift_up(p, p->size, obj);
	p->size++;
	return 0;
}

static void _shift_up(PriorityQueue *p, size_t i, void *o) {
	while (i > 0) {
		size_t parent = (i - 1) >> 1;
		if ((*p->comparator)(p->queue[parent], o) <= 0) {
			break;
		}
		p->queue[i] = p->queue[parent];
		i = parent;
	}
	p->queue[i] = o;
}

static void _shift_down(PriorityQueue *p, size_t i, void *o) {
	size_t h = p->size >> 1;
	while (i < h) {
		size_t ci = (i << 1) + 1;
		size_t r = ci + 1;
		void *c = p->queue[ci];
		if (r < p->size && (*p->comparator)(c, p->queue[r]) > 0) {
			c = p->queue[ci = r];
		}
		if ((*p->comparator)(o, c) <= 0) {
			break;
		}
		p->queue[i] = p->queue[ci];
		i = ci;
	}
	p->queue[i] = o;

}

static int _index_of(PriorityQueue *p, void *o) {
	for (int i = 0; i < p->size; i++) {
		if ((*p->comparator)(p->queue[i], o) == 0) {
			return i;
		}
	}
	return ERR;
}

static int _remove_at(PriorityQueue *p, size_t i) {
	if (i < 0 || i >= p->size) {
		return ERR;
	}
	p->size--;
	if (i == p->size) {
		p->queue[i] = NULL;
	} else {
		void *o = p->queue[i] = p->queue[p->size];
		_shift_down(p, i, o);
		if (o == p->queue[i]) {
			_shift_up(p, i, o);
		}
	}
	_shrink(p);
	return i;
}

static int _shrink(PriorityQueue *p) {
	if (p->capacity <= DEF_SIZE) {
		return ERR;
	}
	if (p->capacity - p->size > p->size) {
		size_t newc = p->size + (p->size >> 1);
		if (newc >= p->capacity) {
			return ERR;
		}
		p->queue = realloc(p->queue, newc * sizeof(void *));
		p->capacity = newc;
	}
	return 0;
}


int remove_element(PriorityQueue *p, void* obj) {
	return _remove_at(p, _index_of(p, obj));
}

void* dequeue(PriorityQueue *p) {
	if (p->size == 0) {
		return NULL;
	}
	void *o = p->queue[0];
	_remove_at(p, 0);
	return o;
}

