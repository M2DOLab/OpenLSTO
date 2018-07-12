// Modification of Zed Shaw's unit testing macros from "Learn C The Hard Way"
// http://c.learncodethehardway.org/book/ex30.html

#undef NDEBUG
#ifndef _MINUNIT_H
#define _MINUNIT_H

#include <stdio.h>
#include <stdlib.h>

#include "debug.h"

/*! \file MinUnit.h
    \brief A simple unit testing framework.
 */

#define mu_suite_start() int message = 0

#define mu_assert(test, message) if (!(test)) { log_err(message); return message; }
#define mu_run_test(test) lsm_debug("\n-----%s", " " #test); \
    message = test(); tests_run++; if (message) return message;

#define RUN_TESTS(name) int main(int argc, char *argv[]) {\
    argc = 1; \
    lsm_debug("----- RUNNING: %s", argv[0]);\
        printf("----\nRUNNING: %s\n", argv[0]);\
        int result = name();\
        if (result != 0) {\
            printf("FAILED: %d\n", result);\
        }\
        else {\
            printf("ALL TESTS PASSED\n");\
        }\
    printf("Tests run: %d\n", tests_run);\
        exit(result != 0);\
}

int tests_run;

#endif /* _MINUNIT_H */
