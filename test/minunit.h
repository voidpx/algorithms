#ifndef TEST_MINUNIT_H_
#define TEST_MINUNIT_H_

#include <stdio.h>
#include <stdlib.h>

/* file: minunit.h */
 #define mu_assert(message, test) do { if (!(test)) { printf("%s\n",message); exit(1);}} while (0)
 #define mu_run_test(test) do { char *message = test(); tests_run++; \
                                if (message) { printf("%s\n",message); exit(1); }} while (0)
 extern int tests_run;


#endif /* TEST_MINUNIT_H_ */
