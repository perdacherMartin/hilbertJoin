
#ifndef MKM_ARGS_H
#define MKM_ARGS_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

void parsing_args(int argc, char* argv[], size_t *n, double *epsilon, size_t *d, size_t *threads, char *filename, bool *isBinary, int *KBLOCK, int *stripes);


#endif //KMEANS_ARGS_H
