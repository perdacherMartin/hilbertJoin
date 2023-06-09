
#ifndef DATA_IO_H
#define DATA_IO_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <omp.h>

#include "mkl_vsl.h"

#define MAX_LINE_LENGTH 2049

int stringEndsWith(const char *str, const char *suffix);
void random_init(double *array, const int N, const int D);
void random_init_8_selective(double *array, const int N, const int D, const int INIT_SEED);
void random_init_unif(double *array, const int N, const int D, const int INIT_SEED);
void read_file(double *array, const int N, const int D, char filename[]);
void save_binary_file(double *array, const int N, const int D, char filename[]);
void save_text_file(double *array, const int N, const int D, char filename[]);

#endif
