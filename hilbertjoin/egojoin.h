#ifndef EGOJOIN_HILBERT_H
#define EGOJOIN_HILBERT_H

#include <omp.h>
#include <immintrin.h>
#include <x86intrin.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "../util/allocation.h"
#include "../util/chrisutil.h"
#include "../util/dataIo.h"
#include "../measure/timer.h"
// #include "../measure/energy.h"
#include "hilloop.h"

#include <boost/lockfree/queue.hpp>
#include <boost/atomic.hpp>

#ifndef KBLOCK
#define KBLOCK 4
#endif

#ifndef NUM_THREADS
#define NUM_THREADS 64
#endif

#define max(a,b) ((a)>(b) ? (a) : (b))
#define min(a,b) ((a)<(b) ? (a) : (b))

struct join_pair {
    size_t p1;
    size_t p2;
};

typedef double vec __attribute__((vector_size(64), aligned(64)));
typedef double vec4 __attribute__((vector_size(32), aligned(32)));
typedef double vecu __attribute__((vector_size(64)));
typedef long long veci64 __attribute__((vector_size(64), aligned(64)));

int epsilonGridCompare(const void *a, const void *b);
void epsilonGridOrdering(size_t n, size_t d, double epsilon, double* array);
int test_ego_loop3(size_t n, size_t d, double epsilon, double *array, long long *result);
int test_ego_loop(size_t n, size_t d, double epsilon, double *array, long long *result);
inline int ceilpowtwo(int m);
void epsilonGridFillList3(size_t n, size_t d, double epsilon, double * array, int *L1, int *L2, int *L3);
void epsilonGridCompleteListMax(size_t n, int *list);
void epsilonGridCompleteListMin(size_t n, int *list);
static inline void transpose_8xd(size_t n, size_t d, double *EGO_array);
void prepareStripes(size_t n, size_t d, int numStripes, double epsilon, double *array, int ** lower, int ** upper, double *self);
static inline long long _custom_mm512_reduce_add_epi64(__m512i a);
static inline void transpose_dx8(size_t n, size_t d, double *EGO_array);
void omp_qsort (void* l, size_t num, size_t size, int (*compar)(const void*,const void*));
// void test_ego_loop3_macro(size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, size_t *countresult, int stripes, int KBLOCK, double *sorttime);
void test_ego_loop3_macro(size_t n, size_t d, double epsilon, double *array, size_t *countresult, int stripes, double *sortTime, double *indextime, double *loadpercent);
// void test_ego_loop3_macro_queue(size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, size_t *countresult, int stripes, int KBLOCK, boost::lockfree::queue<join_pair> &jpartners, double *sorttime);
int cmp_reorder_dim(const void * a, const void *b);
void reorder_dimensions(int n, int d, double *array, int *reorder_dim);
void outputStatistics(int n, int d, double epsilon, double *array, int *reorder_dim);
void sampleHistograms(int n, int d, double epsilon, double *array, int *reorder_dim);

void prepareStripesNoSelf(int nA, int nB, int d, int activeDimensions, double epsilon, double *A, double *B, int ** lower, int **upper, double *selfA, double *selfB);
void test_ego_loop3_noself(const size_t nA, const size_t nB, const int d, const double epsilon, double *arrayA, double *arrayB, size_t *countresult, const int activedims, double *sortTime, double *indextime, double *loadpercent);

#define hilbert_swap(i,j) {memcpy(b, (i), size); memcpy((i), (j), size); memcpy((j), b, size);}
#define hilbert_max(a,b) ((a)>(b) ? (a) : (b))
#define hilbert_min(a,b) ((a)<(b) ? (a) : (b))

extern double EGO_epsilon;
extern int EGO_d;
extern int min_size_qsort;

extern long long * costref;


#define transposeAVX512(r1,r2,r3,r4,r5,r6,r7,r8){\
            vec r1a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r1), mask1, _mm512_castpd_si512 (r2)));\
            vec r2a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r2), mask2, _mm512_castpd_si512 (r1)));\
            vec r3a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r3), mask1, _mm512_castpd_si512 (r4)));\
            vec r4a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r4), mask2, _mm512_castpd_si512 (r3)));\
            vec r5a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r5), mask1, _mm512_castpd_si512 (r6)));\
            vec r6a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r6), mask2, _mm512_castpd_si512 (r5)));\
            vec r7a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r7), mask1, _mm512_castpd_si512 (r8)));\
            vec r8a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r8), mask2, _mm512_castpd_si512 (r7)));\
            vec r1b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r1a), mask3, _mm512_castpd_si512 (r3a)));\
            vec r3b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r3a), mask4, _mm512_castpd_si512 (r1a)));\
            vec r2b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r2a), mask3, _mm512_castpd_si512 (r4a)));\
            vec r4b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r4a), mask4, _mm512_castpd_si512 (r2a)));\
            vec r5b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r5a), mask3, _mm512_castpd_si512 (r7a)));\
            vec r7b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r7a), mask4, _mm512_castpd_si512 (r5a)));\
            vec r6b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r6a), mask3, _mm512_castpd_si512 (r8a)));\
            vec r8b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r8a), mask4, _mm512_castpd_si512 (r6a)));\
            r1 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r1b), mask5, _mm512_castpd_si512 (r5b)));\
            r5 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r5b), mask6, _mm512_castpd_si512 (r1b)));\
            r2 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r2b), mask5, _mm512_castpd_si512 (r6b)));\
            r6 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r6b), mask6, _mm512_castpd_si512 (r2b)));\
            r3 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r3b), mask5, _mm512_castpd_si512 (r7b)));\
            r7 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r7b), mask6, _mm512_castpd_si512 (r3b)));\
            r4 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r4b), mask5, _mm512_castpd_si512 (r8b)));\
            r8 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r8b), mask6, _mm512_castpd_si512 (r4b)));\
            }

#define EGO_PARALLEL(n,d,epsilon,array) {\
    size_t EGO_n = (n);\
    EGO_d = (d);\
    EGO_epsilon = (epsilon);\
    double * EGO_array = (array);\
    int stripes = 2;\
    epsilonGridOrdering(EGO_n, EGO_d, EGO_epsilon, EGO_array);\
    int nn = ceilpowtwo(EGO_n);\
    int **lower = (int **) ddr_alloc (2*sizeof(int*));\
    int **upper = (int **) ddr_alloc (2*sizeof(int*));\
    double *self = (double*) ddr_alloc(sizeof (double) * EGO_n);\
    lower[0] = (int *) ddr_alloc(sizeof (int) * 2 * nn);\
    upper[0] = (int *) ddr_alloc(sizeof (int) * 2 * nn);\
    lower[1] = (int *) ddr_alloc(sizeof (int) * 2 * nn);\
    upper[1] = (int *) ddr_alloc(sizeof (int) * 2 * nn);\
    for (int i=0 ; i<EGO_n ; i++)\
        lower[0][i+nn] = i ;\
    epsilonGridFillList3(EGO_n, EGO_d, EGO_epsilon, EGO_array, upper[0] + nn, lower[1] + nn, upper[1] + nn);\
    for (int j=0 ; j<2 ; j++){\
        for (int i = 0; i < EGO_n; i++)\
           lower[j][nn + i] = lower[j][nn + i] / 8;\
        for (int i = 0; i < EGO_n; i++)\
            upper[j][nn + i] = (upper[j][nn + i] + 7) / 8;\
        epsilonGridCompleteListMin(nn, lower[j]);\
        epsilonGridCompleteListMax(nn, upper[j]);\
    }\
    EGO_epsilon = EGO_epsilon * EGO_epsilon / 2;\
    for (int i=0 ; i<EGO_n ; i++){\
        double h=EGO_epsilon;\
        for(int j=0 ; j<EGO_d ; j++)\
            h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];\
        self[i]=h/2;\
    }\
    double value = -1.;\
    veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};\
    veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};\
    veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};\
    veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};\
    veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};\
    veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};\
    long long overall_load = 0;\
    for (int i = 0; i < EGO_n / 8; i++)\
        for (int j = 0; j < stripes; j++)\
            overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
    int loadstart[NUM_THREADS + 2];\
    for (int i = 0; i <= NUM_THREADS; i++)\
        loadstart[i] = 0;\
    long long cum_load = 0;\
    for (int i = 0; i < EGO_n / 8; i++) {\
        for (int j = 0; j < stripes; j++)\
            cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
        loadstart[cum_load * NUM_THREADS / overall_load + 1] = i;\
    }\
    loadstart[NUM_THREADS] = EGO_n / 8;\
    for (int i = 1; i <= NUM_THREADS; i++)\
        if (loadstart[i] == 0)\
            loadstart[i] = loadstart[i - 1];

#define EGO_PREPARE for (int par = 0; par < NUM_THREADS; par++) {

#define EGO_PARALLEL_TRAN(n,d,epsilon,stripes,array){\
        int EGO_n = (n);\
        int EGO_d = (d);\
        double EGO_epsilon = (epsilon);\
        double * EGO_array = (array);\
        int EGO_blocks = (EGO_d + KBLOCK - 1) / KBLOCK;\
        int EGO_stripes = (stripes);\
        unsigned long long usedload[NUM_THREADS]; for(int i=0; i<NUM_THREADS; i++) usedload[i]=0ull;\
        omp_lock_t criticallock; omp_init_lock(&criticallock); int scritical = -1;\
        CUtilTimer sortTimer;\
        sortTimer.start();\
        epsilonGridOrdering(EGO_n, EGO_d, EGO_epsilon, EGO_array);\
        sortTimer.stop();\
        int nn = ceilpowtwo(EGO_n);\
        int **lower = (int **) malloc (EGO_stripes*sizeof(int*));\
        int **upper = (int **) malloc (EGO_stripes*sizeof(int*));\
        double *self = callocA64(sizeof (double) * EGO_n * EGO_blocks);\
        prepareStripes(EGO_n, EGO_d, EGO_stripes, EGO_epsilon, EGO_array, lower, upper, (double *)0);\
        EGO_epsilon = EGO_epsilon * EGO_epsilon / 2;\
        for (int i=0 ; i<EGO_n ; i++){\
            double h=EGO_epsilon;\
            for(int j=0 ; j<2*KBLOCK && j<EGO_d ; j++)\
                h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];\
            self[i/8*8*EGO_blocks+i%8]=h/2;\
        }\
        for(int k=2 ; k<EGO_blocks ; k++)\
            for (int i=0 ; i<EGO_n ; i++){\
                double h=0;\
                for(int j=k * KBLOCK ; j<(k+1) * KBLOCK && j<EGO_d ; j++)\
                h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];\
            self[(i/8*EGO_blocks+k)*8+i%8]=h/2;\
        }\
        transpose_8xd(EGO_n, EGO_d, EGO_array);\
        unsigned long long *sloadcumcum = (unsigned long long *) malloc(EGO_stripes * (n+7)/8 * sizeof(unsigned long long));\
        sloadcumcum[EGO_stripes-1]=upper[EGO_stripes-1][nn/8]-lower[EGO_stripes-1][nn/8];\
        for(int s=EGO_stripes-1; s>=0 ; s--)\
            sloadcumcum[s] = sloadcumcum[s+1] + upper[s][nn/8] - lower[s][nn/8];\
        for(int i=1; i<(n+7)/8 ; i++)\
            for(int s=EGO_stripes-1; s>=0 ; s--)\
                sloadcumcum[i*EGO_stripes+s] = sloadcumcum[(i-1)*EGO_stripes+s] + upper[s][nn/8+i] - lower[s][nn/8+i];\
        for(int i=0; i<(n+7)/8 ; i++)\
            for(int s=EGO_stripes-2; s>=0 ; s--)\
                sloadcumcum[i*EGO_stripes+s] += sloadcumcum[i*EGO_stripes+s+1];\
        unsigned long long * assignedload = (unsigned long long *) calloc(NUM_THREADS, sizeof (unsigned long long));\
        long long overall_load = 0;\
        for (int i = 0; i < EGO_n / 8; i++)\
            for (int j = 0; j < EGO_stripes; j++)\
                overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
        int *loadstart = (int *) calloc((NUM_THREADS + 1) * EGO_stripes, sizeof(int));\
        long long cum_load = 0;\
        for (int i = 0; i < EGO_n / 8; i++) {\
            for (int j = 0; j < EGO_stripes; j++)\
                cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
            loadstart[cum_load * NUM_THREADS / overall_load + 1] = i;\
        }\
        loadstart[NUM_THREADS] = EGO_n / 8;\
        for (int i = 1; i <= NUM_THREADS; i++)\
            if (loadstart[i] == 0)\
                loadstart[i] = loadstart[i - 1];

#define EGO_LOOP\
        int imin = loadstart[par];\
        int imax = loadstart[par + 1];\
        for (int s = 0; s < stripes; s++) {\
            int i = 0;\
            int j = 0;\
            FGF_HILBERT_FOR(i, j, EGO_n / 8, EGO_n / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],\
                    FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&\
                    FURHIL_ub1 - 1 > lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&\
                    FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {\
                    for (int II = 8 * i; II < 8 * i + 8; II += 2) {\
                        register vec vi1 = _mm512_load_pd(EGO_array + EGO_d * II);\
                        register vec vi2 = _mm512_load_pd(EGO_array + EGO_d * (II + 1));\
                        register vec vj = _mm512_load_pd(EGO_array + EGO_d * 8 * j);\
                        register vec sum1 = _mm512_mul_pd(vi1, vj);\
                        register vec sum9 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 1));\
                        register vec sum2 = _mm512_mul_pd(vi1, vj);\
                        register vec sum10 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 2));\
                        register vec sum3 = _mm512_mul_pd(vi1, vj);\
                        register vec sum11 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 3));\
                        register vec sum4 = _mm512_mul_pd(vi1, vj);\
                        register vec sum12 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 4));\
                        register vec sum5 = _mm512_mul_pd(vi1, vj);\
                        register vec sum13 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 5));\
                        register vec sum6 = _mm512_mul_pd(vi1, vj);\
                        register vec sum14 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 6));\
                        register vec sum7 = _mm512_mul_pd(vi1, vj);\
                        register vec sum15 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 7));\
                        register vec sum8 = _mm512_mul_pd(vi1, vj);\
                        register vec sum16 = _mm512_mul_pd(vi2, vj);\
                        for (register int k = 8; k < d; k += 8) {\
                            vi1 = _mm512_load_pd(EGO_array + EGO_d * II + k);\
                            vi2 = _mm512_load_pd(EGO_array + EGO_d * (II + 1) + k);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * 8 * j + k);\
                            sum1 = _mm512_fmadd_pd(vi1, vj, sum1); \
                            sum9 = _mm512_fmadd_pd(vi2, vj, sum9);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 1) + k);\
                            sum2 = _mm512_fmadd_pd(vi1, vj, sum2);\
                            sum10 = _mm512_fmadd_pd(vi2, vj, sum10);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 2) + k);\
                            sum3 = _mm512_fmadd_pd(vi1, vj, sum3);\
                            sum11 = _mm512_fmadd_pd(vi2, vj, sum11);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 3) + k);\
                            sum4 = _mm512_fmadd_pd(vi1, vj, sum4);\
                            sum12 = _mm512_fmadd_pd(vi2, vj, sum12);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 4) + k);\
                            sum5 = _mm512_fmadd_pd(vi1, vj, sum5);\
                            sum13 = _mm512_fmadd_pd(vi2, vj, sum13);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 5) + k);\
                            sum6 = _mm512_fmadd_pd(vi1, vj, sum6);\
                            sum14 = _mm512_fmadd_pd(vi2, vj, sum14);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 6) + k);\
                            sum7 = _mm512_fmadd_pd(vi1, vj, sum7);\
                            sum15 = _mm512_fmadd_pd(vi2, vj, sum15);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 7) + k);\
                            sum8 = _mm512_fmadd_pd(vi1, vj, sum8);\
                            sum16 = _mm512_fmadd_pd(vi2, vj, sum16);\
                        }\
                        for (int III = 0; III < 2; III++) {\
                            vec indicator;\
                            if (III) {\
                                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);\
                                indicator = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +\
                                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[II + 1]);\
                                if(i==j)\
                                    indicator = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(indicator), __mmask8(255 >> (7-(II&7|1))), *(long long *)&value));\
                            } else {\
                                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);\
                                indicator = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +\
                                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[II]);\
                                if(i==j)\
                                    indicator = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(indicator), __mmask8(255 >> (7-(II&6))), *(long long *)&value));\
                            }

#define EGO_CONSOLIDATE }}} FGF_HILBERT_END(i, j); }

#define EGO_END    }}
#define EGO_END_TRAN    }transpose_dx8(EGO_n, EGO_d, EGO_array);}
#define EGO_END_RESUME  }

#define EGO_LOOP_TRAN\
    int imin = loadstart[par];\
    int imax = loadstart[par + 1];\
    for (int s = 0; s < EGO_stripes; s++) {\
        int i = 0;\
        int j = 0;\
        register veci64 const1 = _mm512_set1_epi64(1LL) ;\
        register veci64 const2 = _mm512_add_epi64(const1, const1) ;\
        register veci64 const3 = _mm512_add_epi64(const1, const2) ;\
        register veci64 const4 = _mm512_add_epi64(const1, const3) ;\
        register veci64 const5 = _mm512_add_epi64(const1, const4) ;\
        register veci64 const6 = _mm512_add_epi64(const1, const5) ;\
        register veci64 const7 = _mm512_add_epi64(const1, const6) ;\
        register veci64 const0 = _mm512_setzero_si512();\
        FGF_HILBERT_FOR(i, j, EGO_n / 8, EGO_n / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],\
                FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&\
                FURHIL_ub1 - 1 >= lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&\
                FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {\
                register vec vi = _mm512_load_pd(self + i * 8 * EGO_blocks);\
                register vec vj = _mm512_load_pd(self + j * 8 * EGO_blocks);\
                register vec sum1 = vi + _mm512_permutexvar_pd(const0, vj);\
                register vec sum2 = vi + _mm512_permutexvar_pd(const1, vj);\
                register vec sum3 = vi + _mm512_permutexvar_pd(const2, vj);\
                register vec sum4 = vi + _mm512_permutexvar_pd(const3, vj);\
                register vec sum5 = vi + _mm512_permutexvar_pd(const4, vj);\
                register vec sum6 = vi + _mm512_permutexvar_pd(const5, vj);\
                register vec sum7 = vi + _mm512_permutexvar_pd(const6, vj);\
                register vec sum8 = vi + _mm512_permutexvar_pd(const7, vj);\
                vi = _mm512_load_pd(EGO_array + (EGO_d * i) * 8);\
                vj = _mm512_load_pd(EGO_array + (EGO_d * j) * 8);\
                sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);\
                sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);\
                sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);\
                sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);\
                sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);\
                sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);\
                sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);\
                sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);\
                int k; for (k=1 ; k<d ; k++){\
                if(k % KBLOCK == 0 && k > KBLOCK){\
                    register veci64 allind = _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63);\
                    allind += _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63);\
                    allind += _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63);\
                    allind += _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63);\
                    allind += _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63);\
                    allind += _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63);\
                    allind += _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63);\
                    allind += _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63);\
                    if(_custom_mm512_reduce_add_epi64(allind) >= 64) {k=d+1; break;}\
                    vi = _mm512_load_pd(self + (i * EGO_blocks + k/KBLOCK) * 8);\
                    vj = _mm512_load_pd(self + (j * EGO_blocks + k/KBLOCK) * 8);\
                    sum1 += vi + _mm512_permutexvar_pd(const0, vj);\
                    sum2 += vi + _mm512_permutexvar_pd(const1, vj);\
                    sum3 += vi + _mm512_permutexvar_pd(const2, vj);\
                    sum4 += vi + _mm512_permutexvar_pd(const3, vj);\
                    sum5 += vi + _mm512_permutexvar_pd(const4, vj);\
                    sum6 += vi + _mm512_permutexvar_pd(const5, vj);\
                    sum7 += vi + _mm512_permutexvar_pd(const6, vj);\
                    sum8 += vi + _mm512_permutexvar_pd(const7, vj);\
                }\
                vi = _mm512_load_pd(EGO_array + (EGO_d * i + k) * 8);\
                vj = _mm512_load_pd(EGO_array + (EGO_d * j + k) * 8);\
                sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);\
                sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);\
                sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);\
                sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);\
                sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);\
                sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);\
                sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);\
                sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);\
            }\
            if(k<=d){{\
                if (i==j){\
                    sum1 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum1), 255, 0xbff0000000000000ull));\
                    sum2 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum2), 254, 0xbff0000000000000ull));\
                    sum3 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum3), 252, 0xbff0000000000000ull));\
                    sum4 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum4), 248, 0xbff0000000000000ull));\
                    sum5 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum5), 240, 0xbff0000000000000ull));\
                    sum6 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum6), 224, 0xbff0000000000000ull));\
                    sum7 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum7), 192, 0xbff0000000000000ull));\
                    sum8 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum8), 128, 0xbff0000000000000ull));\
                }


#define EGO_PARALLEL_TRAN_NOSELF(nA,nB,d,epsilon,activedim,arrayA,arrayB) {\
    int EGO_nA = (nA);\
    int EGO_nB = (nB);\
    int EGO_d = (d);\
    double EGO_epsilon = (epsilon);\
    double * EGO_arrayA = (arrayA);\
    double * EGO_arrayB = (arrayB);\
    int EGO_blocks = (EGO_d + KBLOCK - 1) / KBLOCK;\
    int EGO_activedim = (activedim);\
    int EGO_stripes = 1; for(int i=0; i<EGO_activedim; i++) EGO_stripes*=3;\
    unsigned long long usedload[NUM_THREADS]; for(int i=0; i<NUM_THREADS; i++) usedload[i]=0ull;\
    omp_lock_t criticallock; omp_init_lock(&criticallock); int scritical = -1;\
    CUtilTimer sortTimer;\
    sortTimer.start();\
    epsilonGridOrdering(EGO_nA, EGO_d, EGO_epsilon, EGO_arrayA);\
    epsilonGridOrdering(EGO_nB, EGO_d, EGO_epsilon, EGO_arrayB);\
    sortTimer.stop();\
    int nn = ceilpowtwo(EGO_nA);\
    int **lower = (int **) malloc (EGO_stripes*sizeof(int*));\
    int **upper = (int **) malloc (EGO_stripes*sizeof(int*));\
    double *selfA = callocA64(sizeof (double) * EGO_nA * EGO_blocks);\
    double *selfB = callocA64(sizeof (double) * EGO_nB * EGO_blocks);\
    prepareStripesNoSelf(EGO_nA, EGO_nB, EGO_d, EGO_activedim, EGO_epsilon, EGO_arrayA, EGO_arrayB, lower, upper, (double *)0, (double*)0);\
    EGO_epsilon = EGO_epsilon * EGO_epsilon / 2;\
    for(int i=EGO_nA ; i<(EGO_nA+7)/8*8 ; i++)\
        for(int j=0 ; j<EGO_d ; j++)\
            EGO_arrayA[i*EGO_d+j] = 1e150 - 1e140*(double)(i%8);\
    for(int i=EGO_nB ; i<(EGO_nB+7)/8*8 ; i++)\
        for(int j=0 ; j<EGO_d ; j++)\
            EGO_arrayB[i*EGO_d+j] = 1e150 - 1e140*(double)(i%8);\
    for (int i=0 ; i<EGO_nA ; i++){\
        double h=EGO_epsilon;\
        for(int j=0 ; j<2*KBLOCK && j<EGO_d ; j++)\
            h-=EGO_arrayA[i*EGO_d+j]*EGO_arrayA[i*EGO_d+j];\
        selfA[i/8*8*EGO_blocks+i%8]=h/2;\
    }\
    for (int i=0 ; i<EGO_nB ; i++){\
        double h=EGO_epsilon;\
        for(int j=0 ; j<2*KBLOCK && j<EGO_d ; j++)\
            h-=EGO_arrayB[i*EGO_d+j]*EGO_arrayB[i*EGO_d+j];\
        selfB[i/8*8*EGO_blocks+i%8]=h/2;\
    }\
    for(int k=2 ; k<EGO_blocks ; k++)\
        for (int i=0 ; i<EGO_nA ; i++){\
            double h=0;\
            for(int j=k * KBLOCK ; j<(k+1) * KBLOCK && j<EGO_d ; j++)\
                h-=EGO_arrayA[i*EGO_d+j]*EGO_arrayA[i*EGO_d+j];\
        selfA[(i/8*EGO_blocks+k)*8+i%8]=h/2;\
    }\
    for(int k=2 ; k<EGO_blocks ; k++)\
        for (int i=0 ; i<EGO_nB ; i++){\
            double h=0;\
            for(int j=k * KBLOCK ; j<(k+1) * KBLOCK && j<EGO_d ; j++)\
                h-=EGO_arrayB[i*EGO_d+j]*EGO_arrayB[i*EGO_d+j];\
        selfB[(i/8*EGO_blocks+k)*8+i%8]=h/2;\
    }\
    transpose_8xd(EGO_nA, EGO_d, EGO_arrayA);\
    if(EGO_arrayA != EGO_arrayB)\
        transpose_8xd(EGO_nB, EGO_d, EGO_arrayB);\
    unsigned long long *sloadcumcum = (unsigned long long *) malloc(EGO_stripes * (EGO_nA+7)/8 * sizeof(unsigned long long));\
    sloadcumcum[EGO_stripes-1]=upper[EGO_stripes-1][nn/8]-lower[EGO_stripes-1][nn/8];\
    for(int s=EGO_stripes-1; s>=0 ; s--)\
        sloadcumcum[s] = sloadcumcum[s+1] + upper[s][nn/8] - lower[s][nn/8];\
    for(int i=1; i<(EGO_nA+7)/8 ; i++)\
        for(int s=EGO_stripes-1; s>=0 ; s--)\
            sloadcumcum[i*EGO_stripes+s] = sloadcumcum[(i-1)*EGO_stripes+s] + upper[s][nn/8+i] - lower[s][nn/8+i];\
    for(int i=0; i<(EGO_nA+7)/8 ; i++)\
        for(int s=EGO_stripes-2; s>=0 ; s--)\
            sloadcumcum[i*EGO_stripes+s] += sloadcumcum[i*EGO_stripes+s+1];\
    unsigned long long * assignedload = (unsigned long long *) calloc(NUM_THREADS, sizeof (unsigned long long));\
    long long overall_load = 0;\
    for (int i = 0; i < EGO_nA / 8; i++)\
        for (int j = 0; j < EGO_stripes; j++)\
            overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
    int *loadstart = (int *) calloc((NUM_THREADS + 1) * EGO_stripes, sizeof(int));\
    long long cum_load = 0;\
    for (int i = 0; i < EGO_nA / 8; i++) {\
        for (int j = 0; j < EGO_stripes; j++)\
            cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
        loadstart[cum_load * NUM_THREADS / overall_load + 1] = i;\
    }\
    loadstart[NUM_THREADS] = EGO_nA / 8;\
    for (int i = 1; i <= NUM_THREADS; i++)\
        if (loadstart[i] == 0)\
            loadstart[i] = loadstart[i - 1];


#define EGO_LOOP_TRAN_NOSELF\
        int imin = loadstart[par];\
        int imax = loadstart[par + 1];\
        for (int s = 0; s < EGO_stripes; s++) {\
            int i = 0;\
            int j = 0;\
            register veci64 const1 = _mm512_set1_epi64(1LL) ;\
            register veci64 const2 = _mm512_add_epi64(const1, const1) ;\
            register veci64 const3 = _mm512_add_epi64(const1, const2) ;\
            register veci64 const4 = _mm512_add_epi64(const1, const3) ;\
            register veci64 const5 = _mm512_add_epi64(const1, const4) ;\
            register veci64 const6 = _mm512_add_epi64(const1, const5) ;\
            register veci64 const7 = _mm512_add_epi64(const1, const6) ;\
            register veci64 const0 = _mm512_setzero_si512();\
            FGF_HILBERT_FOR(i, j, EGO_nA / 8, EGO_nB / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],\
                FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&\
                FURHIL_ub1 - 1 >= lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&\
                FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {\
                    register vec vi = _mm512_load_pd(selfA + i * 8 * EGO_blocks);\
                    register vec vj = _mm512_load_pd(selfB + j * 8 * EGO_blocks);\
                    register vec sum1 = vi + _mm512_permutexvar_pd(const0, vj);\
                    register vec sum2 = vi + _mm512_permutexvar_pd(const1, vj);\
                    register vec sum3 = vi + _mm512_permutexvar_pd(const2, vj);\
                    register vec sum4 = vi + _mm512_permutexvar_pd(const3, vj);\
                    register vec sum5 = vi + _mm512_permutexvar_pd(const4, vj);\
                    register vec sum6 = vi + _mm512_permutexvar_pd(const5, vj);\
                    register vec sum7 = vi + _mm512_permutexvar_pd(const6, vj);\
                    register vec sum8 = vi + _mm512_permutexvar_pd(const7, vj);\
                    vi = _mm512_load_pd(EGO_arrayA + (EGO_d * i) * 8);\
                    vj = _mm512_load_pd(EGO_arrayB + (EGO_d * j) * 8);\
                    sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);\
                    sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);\
                    sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);\
                    sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);\
                    sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);\
                    sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);\
                    sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);\
                    sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);\
                    int k; for (k=1 ; k<d ; k++){\
                        if(k % KBLOCK == 0 && k > KBLOCK){\
                            register veci64 allind = _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63);\
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63);\
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63);\
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63);\
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63);\
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63);\
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63);\
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63);\
                            if(_custom_mm512_reduce_add_epi64(allind) >= 64) {k=d+1; break;}\
                            vi = _mm512_load_pd(selfA + (i * EGO_blocks + k/KBLOCK) * 8);\
                            vj = _mm512_load_pd(selfB + (j * EGO_blocks + k/KBLOCK) * 8);\
                            sum1 += vi + _mm512_permutexvar_pd(const0, vj);\
                            sum2 += vi + _mm512_permutexvar_pd(const1, vj);\
                            sum3 += vi + _mm512_permutexvar_pd(const2, vj);\
                            sum4 += vi + _mm512_permutexvar_pd(const3, vj);\
                            sum5 += vi + _mm512_permutexvar_pd(const4, vj);\
                            sum6 += vi + _mm512_permutexvar_pd(const5, vj);\
                            sum7 += vi + _mm512_permutexvar_pd(const6, vj);\
                            sum8 += vi + _mm512_permutexvar_pd(const7, vj);\
                       }\
                       vi = _mm512_load_pd(EGO_arrayA + (EGO_d * i + k) * 8);\
                       vj = _mm512_load_pd(EGO_arrayB + (EGO_d * j + k) * 8);\
                       sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);\
                       sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);\
                       sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);\
                       sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);\
                       sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);\
                       sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);\
                       sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);\
                       sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);\
                   }\
                   if(k<=d){{

#define EGO_END_TRAN_NOSELF    }transpose_dx8(EGO_nA, EGO_d, EGO_arrayA);if(EGO_arrayA!=EGO_arrayB)transpose_dx8(EGO_nB, EGO_d, EGO_arrayB);}



#endif
