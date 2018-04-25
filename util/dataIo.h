
#ifndef DATA_IO_H
#define DATA_IO_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <omp.h>

#define MAX_LINE_LENGTH 2049

void random_init(double *array, const int N, const int D);
void read_file(double *array, const int N, const int D, char filename[], const bool IS_BINARY);
void save_binary_file(double *array, const int N, const int D, char filename[]);
void save_text_file(double *array, const int N, const int D, char filename[]);

#endif
