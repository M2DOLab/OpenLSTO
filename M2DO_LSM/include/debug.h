// Zed Shaw's debug macros taken from "Learn C The Hard Way"
// http://c.learncodethehardway.org/book/ex20.html

#ifndef _DEBUG_H
#define _DEBUG_H

#include <stdio.h>
#include <errno.h>
#include <string.h>

/*! \file Debug.h
    \brief A set of debugging macros.
 */

#ifdef NDEBUG
#define lsm_debug(M, ...)
#else
#define lsm_debug(M, ...) fprintf(stderr, "[DEBUG] %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif

#define lsm_clean_errno() (errno == 0 ? "None" : strerror(errno))

#define lsm_log_err(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, lsm_clean_errno(), ##__VA_ARGS__)

#define lsm_log_warn(M, ...) fprintf(stderr, "[WARN] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, lsm_clean_errno(), ##__VA_ARGS__)

#define lsm_log_info(M, ...) fprintf(stderr, "[INFO] (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#define lsm_check(A, M, ...) if(!(A)) { lsm_log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define lsm_sentinel(M, ...)  { lsm_log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define lsm_check_mem(A) check((A), "Out of memory.")

#define lsm_check_debug(A, M, ...) if(!(A)) { lsm_debug(M, ##__VA_ARGS__); errno=0; goto error; }

#endif /* _DEBUG_H */
