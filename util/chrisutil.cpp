
#include "chrisutil.h"

double * mallocA64(size_t s) {
    long long adr = (long long) malloc(s + 72);
    long long adr2 = (adr + 71) & ~63;
    ((long long *) adr2)[-1] = adr;
    return (double *) adr2;
}

double * callocA64(size_t s) {
    long long adr = (long long) calloc(s + 72, 1);
    long long adr2 = (adr + 71) & ~63;
    ((long long *) adr2)[-1] = adr;
    return (double *) adr2;
}

void freeA64(void * adr) {
    free((void *) (((long long *) adr)[-1]));
}

double stopc(void) {
    static double last;
    double c = clock();
    double r = c - last;
    last = c;
    return r / CLOCKS_PER_SEC;
}

double timestamp(void){
    struct timeval c;
    gettimeofday(&c, NULL);
    long long r = (c.tv_sec * 100 + c.tv_usec / 10000) % 1000000000000000 ;
    return (double)r / 100.;
}

double stop(void) {
    static struct timeval last;
    struct timeval c;
    gettimeofday(&c, NULL);
    double r = (double) (c.tv_sec - last.tv_sec) + (double) (c.tv_usec - last.tv_usec) / 1000000.;
    last.tv_sec = c.tv_sec;
    last.tv_usec = c.tv_usec;
    return r;
}
