#ifndef ALLOCATION_H
#define ALLOCATION_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

#ifdef __APPLE__
  #include <cstdlib>
  #include <malloc/malloc.h>
#else
  #if defined(__INTEL_COMPILER)
    #include <malloc.h>
  #else
    #include <mm_malloc.h>
  #endif
#endif

#include <hbwmalloc.h>
#include <errno.h>

#define ALIGNMENT 64

void * ddr_alloc(size_t bytes);
void ddr_free(void * ptrs);

void ** ompx_ddr_calloc(size_t bytes);
void ompx_ddr_free(void ** ptrs);

void ompx_hbm_errcode_check(int errcode);
void ** ompx_hbm_calloc(size_t bytes);
void ompx_hbm_free(void ** ptrs);


#endif
