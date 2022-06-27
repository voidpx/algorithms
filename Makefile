
RM := rm -rf

CFLAGS = -O0 -g

all: pq.so crc test

crc: src/crc.c
	@echo 'building $@...'
	gcc $(CFLAGS) -o $@ $^

pq.so: src/priorityqueue.c
	@echo 'building $@...'
	gcc $(CFLAGS) -fPIC -shared -o $@ $^

clean:
	@echo 'cleaning...'
	-$(RM) pq.so crc
	-@echo ' '

test:
	@echo 'testing...'
	gcc -o test_pq test/test_pq.c src/priorityqueue.c
	chmod +x ./test_pq && ./test_pq
	$(RM) ./test_pq
	@echo 'testing priority queue successful'

.PHONY: all clean test
