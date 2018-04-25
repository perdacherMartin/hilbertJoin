#include "allocation.h"

void *ddr_alloc(size_t bytes){
    void *ptr=NULL;
#ifdef __APPLE__
    if ( posix_memalign((void **)&ptr, ALIGNMENT, bytes) != 0 ) {
        fprintf(stderr, "Error in allocating memory with ddr_alloc!\n");
        exit(1);
    }
#else
    ptr = _mm_malloc(bytes, ALIGNMENT);

    if ( ptr == NULL ){
        fprintf(stderr, "Error in allocating memory with ddr_alloc!\n");
        exit(1);
    }

#endif

    return ptr;
}

void ddr_free(void *ptrs){

#ifdef __APPLE__
    free(ptrs);
#else
    _mm_free(ptrs);
#endif

}

// quadrant mode, nested omp, ddr allocation
void ** ompx_ddr_calloc(size_t bytes){

    int np = omp_get_max_threads(); // returns for nested 4,64 --> 4
    void ** ptrs = (void**) _mm_malloc(np * sizeof(void*),ALIGNMENT);

    // printf("omp_get_max_threads: %d\n", np);
    if ( ptrs == NULL ){
        fprintf(stderr, "Error in allocating ddr memory!\n");
        exit(1);
    }

    #pragma omp parallel shared(ptrs)
    {
        int me = omp_get_thread_num();
        ptrs[me] = _mm_malloc((bytes / np) + 1,ALIGNMENT);
        memset(ptrs[me], 0, (bytes / np) + 1);
    }

    return ptrs;
}

// quadrant mode, nested omp, ddr free
void ompx_ddr_free(void ** ptrs){
    int np = omp_get_max_threads();
    printf("omp_get_max_threads: %d, \n", np);
    #pragma omp parallel shared(ptrs)
    {
        int me = omp_get_thread_num();
        _mm_free(ptrs[me]);
    }
    _mm_free(ptrs);
}



// quadrant mode, hbm allocation error handling
void ompx_hbm_errcode_check(int errcode){
    if ( errcode != 0 ){
        if ( errcode == EINVAL ){
            fprintf(stderr, "Error in allocating hbm memory!\n Alginment parameter is not a multiple of 2.\n");
            exit(1);
        }

        if ( errcode == ENOMEM ){
            fprintf(stderr, "Error in allocating hbm memory!\n There is insufficient memory to satisfy the request.\n");
            exit(1);
        }
    }
}

// quadrant mode, nested omp, hbm calloc
void ** ompx_hbm_calloc(size_t bytes){
    int np = omp_get_max_threads(); // returns for nested 4,64 --> 4
    void ** ptrs = NULL;

    int errcode = hbw_posix_memalign((void**) &ptrs,ALIGNMENT,np * sizeof(void*));

    ompx_hbm_errcode_check(errcode);

    #pragma omp parallel shared(ptrs)
    {
        int me = omp_get_thread_num();
        int all = omp_get_max_threads();
        errcode = hbw_posix_memalign((void**) &ptrs[me],ALIGNMENT,(bytes/all) + 1);
        ompx_hbm_errcode_check(errcode);
        memset(ptrs[me], 0, (bytes/all) + 1);
    }

    return ptrs;
}

// nested omp, hbm free
void ompx_hbm_free(void ** ptrs){
    int np = omp_get_max_threads();
    // printf("omp_get_max_threads: %d, \n", np);
    #pragma omp parallel shared(ptrs)
    {
        int me = omp_get_thread_num();
        hbw_free(ptrs[me]);
    }
    hbw_free(ptrs);
}
