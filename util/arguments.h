
#ifndef MKM_ARGS_H
#define MKM_ARGS_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

void parsing_args(int argc, char* argv[], size_t *n, double *epsilon, size_t *d, char *filename, int *stripes);
void parsing_args_join(int argc, char* argv[], size_t *n, size_t *m, double *epsilon, size_t *d, char *filename, char *filename2, int *stripes);


#endif //KMEANS_ARGS_H
