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
