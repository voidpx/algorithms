#ifndef PRIORITYQUEUE_H_
#define PRIORITYQUEUE_H_

#include <stddef.h>

#define ERR -1;

typedef struct _PriorityQueue {
	void **queue;
	size_t capacity;
	size_t size;
	int (*comparator)(void *, void *);
	void *lock; // could be used for synchronization, unused for now
} PriorityQueue;

PriorityQueue* new_pq(size_t capacity, int (*comparator)(void *, void *));
void free_pq(PriorityQueue *p);
int enqueue(PriorityQueue *p, void *obj);
int remove_element(PriorityQueue *p, void *obj);
void* dequeue(PriorityQueue *p);

#endif /* PRIORITYQUEUE_H_ */
