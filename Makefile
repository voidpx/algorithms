
RM := rm -rf
CC = gcc
CFLAGS = -O0 -g

all: pq.so crc fft.so test

crc: src/crc.c
	@echo 'building $@...'
	$(CC) $(CFLAGS) -o $@ $^

pq.so: src/priorityqueue.c
	@echo 'building $@...'
	$(CC) $(CFLAGS) -fPIC -shared -o $@ $^

fft.so: src/fft.c
	@echo 'building $@...'
	$(CC) $(CFLAGS) -fPIC -shared -o $@ $^ -lm

clean:
	@echo 'cleaning...'
	-$(RM) fft.so pq.so crc
	-@echo ' '
 
TEST_PQ = test/test_pq
TEST_FFT = test/test_fft
 
test:
	@echo 'testing...'
	$(CC) $(CFLAGS) -o $(TEST_PQ) test/test_pq.c src/priorityqueue.c
	chmod +x $(TEST_PQ) && $(TEST_PQ)
	$(RM) $(TEST_PQ)
	@echo 'testing priority queue successful'
	
	@echo 'testing fft'
	$(CC) $(CFLAGS) -g -o $(TEST_FFT) test/test_fft.c src/fft.c -lm
	chmod +x $(TEST_FFT) && $(TEST_FFT)
	$(RM) $(TEST_FFT)
	@echo 'testing fft successful'
	
.PHONY: all clean test
