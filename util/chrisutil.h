#ifndef CHRIS_UTIL_H
#define CHRIS_UTIL_H

#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#define DBL_MAX         1.7976931348623158e+308

double timestamp(void);
double * mallocA64(size_t s);
double * callocA64(size_t s);
void freeA64(void * adr);
double stopc(void);
double stop(void);

#endif
