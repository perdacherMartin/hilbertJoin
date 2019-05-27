
#include "egojoin.h"

double EGO_epsilon;
int EGO_d;
int min_size_qsort;
long long * costref;

int epsilonGridCompare(const void *a, const void *b) {
    double * A = (double *) a;
    double * B = (double *) b;
    for (int i = 0; i < EGO_d; i++) {
        int h = (int) (floor(A[i] / EGO_epsilon) - floor(B[i] / EGO_epsilon));
        if (h != 0) return h;
    }
    return 0;
}

void sampleHistograms(int n, int d, double epsilon, double *array, int *reorder_dim) {
    double mins[d];
    double maxs[d];
    memcpy(mins, array, d * sizeof (double));
    memcpy(maxs, array, d * sizeof (double));
#pragma omp parallel for reduction(min:mins[:d]) reduction(max:maxs[:d])
    for (int par = 0; par < NUM_THREADS; par++) {
        int i = par * n / NUM_THREADS;
        memcpy(mins, array + i*d, d * sizeof (double));
        memcpy(maxs, array + i*d, d * sizeof (double));
//        for (i++; i < (par + 1) * n / NUM_THREADS; i += (int) (-log(1. - drand48())*100.)) {
        for (i++; i < (par + 1) * n / NUM_THREADS; i++) {
            for (int j = 0; j < d; j++) {
                mins[j] = min(mins[j], array[i * d + j]);
                maxs[j] = max(maxs[j], array[i * d + j]);
            }
        }
    }
    int sizes[d];
    int cum = 0;
    int offs[d + 1];
    for (int j = 0; j < d; j++) {
        mins[j] = floor(mins[j] / epsilon);
        maxs[j] = floor(maxs[j] / epsilon);
        sizes[j] = (int) maxs[j] - (int) mins[j] + 1;
        offs[j] = cum - (int) mins[j];
        cum += sizes[j];
    }
    offs[d] = cum;
    long long * histos = (long long *) calloc(cum, sizeof (long long));
#pragma omp parallel for reduction(+:histos[:cum])
    for (int par = 0; par < NUM_THREADS; par++)
//        for (int i = par * n / NUM_THREADS; i < (par + 1) * n / NUM_THREADS; i += (int) (-log(1. - drand48())*100.))
        for (int i = par * n / NUM_THREADS; i < (par + 1) * n / NUM_THREADS; i++)
            for (int j = 0; j < d; j++)
                histos[(int) min(maxs[j],max(mins[j], floor(array[i * d + j] / epsilon))) + offs[j]]++;
    int h = 0;
    costref = (long long *) calloc(d, sizeof (long long));
    for (int i = 0; i < d; i++) {
#ifdef OUTPUT
        printf("%2d ", i);
#endif
        for (int j = 0; j < sizes[i]; j++) {
            costref[i] += histos[h]*(histos[h] - 1) / 2 + (j > 0 ? histos[h - 1] * histos[h] : 0);
#ifdef OUTPUT
            printf("%lld ", histos[h]);
#endif
            h++;
        }
#ifdef OUTPUT
        printf(" => %lld\n", costref[i]);
#endif
    }
    int * reorder_rev = (int*) malloc((d + 8) * sizeof (int));
    for (int j = 0; j < d + 8; j++)
        reorder_rev[j] = j;
    qsort(reorder_rev, d, sizeof (int), cmp_reorder_dim);
#ifdef OUTPUT
    for (int j = 0; j < d + 8; j++)
        printf("%2d %2d %lld\n", j, reorder_rev[j], j < d ? costref[reorder_rev[j]] : 0);
#endif
    // reorder_dim = (int*) malloc((d + 8) * sizeof (int));
    for (int j = 0; j < d + 8; j++)
        reorder_dim[reorder_rev[j]] = j;
    free(reorder_rev);
    free(costref);
    free(histos);
}

void test_ego_loop(size_t n, size_t d, double epsilon, double *array, long long *result){
    // int KBLOCK=8;
    *result = 0l;
    long long iresult = 0l;
    EGO_PARALLEL(n, d, epsilon, array)
        #pragma omp parallel for reduction(+: iresult)
    EGO_PREPARE
        // printf("parallel %d\n", omp_get_thread_num() );
        veci64 resultvec = _mm512_setzero_si512();
    EGO_LOOP{
        resultvec += _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(indicator), 63));
    }
    EGO_CONSOLIDATE{
        long long resultvecX[8];
        _mm512_storeu_si512((void *) resultvecX, resultvec);
        //    printf("%ld %ld %ld %ld %ld %ld %ld %ld\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
        //    printf("%lx %lx %lx %lx %lx %lx %lx %lx\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
        iresult += resultvecX[0] + resultvecX[1] + resultvecX[2] + resultvecX[3] + resultvecX[4] + resultvecX[5] + resultvecX[6] + resultvecX[7];
//        printf("par = %d: %d %d\n", par, result, _mm512_reduce_add_epi64(resultvec));
      }
    EGO_END

#ifdef OUTPUT
    printf("count of join partners: %lld\n", iresult);
#endif
    *result = iresult;
}

void test_ego_loop3(size_t n, size_t d, double epsilon, double *array, long long *result){
    *result = 0l;
    // int KBLOCK=8;
    long long iresult=0l;
    omp_set_num_threads(NUM_THREADS);
    EGO_PARALLEL_TRAN(n, d, epsilon, 5, array)
        #pragma omp parallel for reduction(+: iresult)
    EGO_PREPARE
        veci64 resultvec = _mm512_setzero_si512();
        veci64 eights = _mm512_set1_epi64(8ll) ;
    EGO_LOOP_TRAN{
        resultvec += eights - _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63) ;
    }
    EGO_CONSOLIDATE{
        iresult += _custom_mm512_reduce_add_epi64(resultvec);
//        double testres[8] __attribute__((aligned(64)));
//        _mm512_store_epi64(testres, resultvec);
//        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
    }
    EGO_END_TRAN
#ifdef OUTPUT
    printf("result %lld\n", iresult);
#endif
    *result = iresult;
}

inline int ceilpowtwo(int m) {
    // determines the next power of two which is >=m, e.g. ceilpowtwo(13) = 16
    m--;
    m |= m >> 1;
    m |= m >> 2;
    m |= m >> 4;
    m |= m >> 8;
    return (m | (m >> 16)) + 1;
}

void epsilonGridOrdering(size_t n, size_t d, double epsilon, double* array) {
    EGO_epsilon = epsilon;
    EGO_d = d;
    min_size_qsort = n/NUM_THREADS*2 ;
#pragma omp parallel
#pragma omp master
    omp_qsort(array, n, d * sizeof (double), epsilonGridCompare);
#pragma omp taskwait

}


void epsilonGridFillList3(size_t n, size_t d, double epsilon, double * array, int *L1, int *L2, int *L3) {
    int a1, a2, a3;
    a1 = a2 = a3 = 0;
    int back1 = 0, back2 = 0, back3 = 0;
    EGO_epsilon = epsilon;
    EGO_d = d;
    double h1[d], h2[d], h3[d];
    for (int b = 0; b < n; b++) {
        for (int i = 0; i < d; i++)
            h1[i] = h3[i] = array[d * b + i] + epsilon;
        h1[0] = array[d * b];
        for (int i = 1; i < d; i++)
            h2[i] = array[d * b + i] - epsilon;
        h2[0] = array[d * b] + epsilon;
        while (a1 > 0 && epsilonGridCompare(array + a1 * d, h1) > 0) {
            a1--;
            back1++;
        }
        while (a1 < n && epsilonGridCompare(array + a1 * d, h1) <= 0)
            a1++;
        while (a2 > 0 && epsilonGridCompare(array + a2 * d, h2) >= 0) {
            a2--;
            back2++;
        }
        while (a2 < n && epsilonGridCompare(array + a2 * d, h2) < 0)
            a2++;
        while (a3 > 0 && epsilonGridCompare(array + a3 * d, h3) > 0) {
            a3--;
            back3++;
        }
        while (a3 < n && epsilonGridCompare(array + a3 * d, h3) <= 0)
            a3++;
        L1[b] = a1;
        L2[b] = a2;
        L3[b] = a3;
    }

#ifdef OUTPUT
    printf("Backsteps: %d %d %d\n", back1, back2, back3);
#endif

}

void epsilonGridCompleteListMax(size_t n, int *list) {
    int m = ceilpowtwo(n);
    for (size_t i = n; i < m; i++)
        list[i+m] = (n+7)/8;
    if (n % 2)
        list[n] = list[2 * n - 1];
    for (size_t i = n - 1; i > 0; i--)
        list[i] = hilbert_max(list[2 * i], list[2 * i + 1]);
}

void epsilonGridCompleteListMin(size_t n, int *list) {
    int m = ceilpowtwo(n);
    for (size_t i = n; i < m; i++)
        list[i+m] = 0;
    if (n % 2)
        list[n] = list[2 * n - 1];
    for (size_t i = n - 1; i > 0; i--)
        list[i] = hilbert_min(list[2 * i], list[2 * i + 1]);
}

static inline void transpose_8xd(size_t n, size_t d, double *EGO_array) {
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++) {
        veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};
        veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};
        veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};
        veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};
        veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};
        veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};
        double help[(d+7) * 8]__attribute__((aligned(64)));
        for (int i = par * n / NUM_THREADS / 8; i < (par + 1) * n / NUM_THREADS / 8; i++) {
            for (int k = 0; k < d; k += 8) {
                vec v1 = _mm512_loadu_pd(EGO_array + 8*i * d + k);
                vec v2 = _mm512_loadu_pd(EGO_array + (8*i + 1) * d + k);
                vec v3 = _mm512_loadu_pd(EGO_array + (8*i + 2) * d + k);
                vec v4 = _mm512_loadu_pd(EGO_array + (8*i + 3) * d + k);
                vec v5 = _mm512_loadu_pd(EGO_array + (8*i + 4) * d + k);
                vec v6 = _mm512_loadu_pd(EGO_array + (8*i + 5) * d + k);
                vec v7 = _mm512_loadu_pd(EGO_array + (8*i + 6) * d + k);
                vec v8 = _mm512_loadu_pd(EGO_array + (8*i + 7) * d + k);
                transposeAVX512(v1, v2, v3, v4, v5, v6, v7, v8);
                _mm512_store_pd(help + k * 8, v1);
                _mm512_store_pd(help + k * 8 + 8, v2);
                _mm512_store_pd(help + k * 8 + 16, v3);
                _mm512_store_pd(help + k * 8 + 24, v4);
                _mm512_store_pd(help + k * 8 + 32, v5);
                _mm512_store_pd(help + k * 8 + 40, v6);
                _mm512_store_pd(help + k * 8 + 48, v7);
                _mm512_store_pd(help + k * 8 + 56, v8);
            }
            memcpy(EGO_array + i * 8 * d, help, d * 8 * sizeof (double));
        }
    }
}


void prepareStripes(size_t n, size_t d, int numStripes, double epsilon, double *array, int ** lower, int ** upper, double *self){
//                    for (int a=0 ; a<d && a<20 ; a++)
//                        printf("%5.3lf ", array[17*1024*d+a]);
//                    printf("\n");
//                    for (int a=0 ; a<d && a<d ; a++)
//                        printf("%3.1lf ", (double)(int)(array[17*1024*d+a]/epsilon)*epsilon);
//                    printf("\n");
    double lowerDisplace[numStripes][d];
    double upperDisplace[numStripes][d];
    int neun = 1;
    int activeDimensions = 0;
    while (neun < 2*numStripes){
        neun*=3 ;
        activeDimensions++;
    }
//    printf("activeDimensions: %d, %d\n", activeDimensions, neun);
    for (int t=0 ; t<numStripes-1 ; t++){
        int k=numStripes+t+1;
        int i=(k*neun*2+2*numStripes)/numStripes/4;
        int i2=((k-1)*neun*2+2*numStripes)/numStripes/4;
        if ((i + 3) / 3 > (i2 + 3) / 3)
            i -= i % 3;
        int umrechnung = 0;
        int power = 1;
        for (int j=i; j>0 ; j/=3){
            umrechnung += power*(j%3) ;
            power *=10;
        }
        // printf("%d %d %d\n", t, i, umrechnung);
        i2 = i;
        for (int j = activeDimensions - 1; j >= 0; j--, i2 /= 3) {
            lowerDisplace[t + 1][j] = epsilon * (i2 % 3 - 1);
        }
        for (int j = activeDimensions; j < d; j++)
            lowerDisplace[t + 1][j] = (-epsilon);
        i2 = i - 1;
        for (int j = activeDimensions-1; j >= 0; j--, i2 /= 3) {
            upperDisplace[t][j] = epsilon * (i2 % 3 - 1);
        }
        for (int j = activeDimensions; j < d; j++)
            upperDisplace[t][j] = epsilon;
    }
    for (int j = 0; j < d; j++)
        lowerDisplace[0][j] = 0.0;
    for(int j=0; j<d ; j++)
        upperDisplace[numStripes-1][j] = epsilon;
//    for(int i=0 ; i<numStripes ; i++){
//        printf("%4.1lf %4.1lf %4.1lf %4.1lf %4.1lf     %4.1lf %4.1lf %4.1lf %4.1lf %4.1lf\n",
//                lowerDisplace[i][0],lowerDisplace[i][1],lowerDisplace[i][2],lowerDisplace[i][3],lowerDisplace[i][4],
//                upperDisplace[i][0],upperDisplace[i][1],upperDisplace[i][2],upperDisplace[i][3],upperDisplace[i][4]) ;
//    }
    int nn = ceilpowtwo(n);
    lower[0] = (int *) callocA64(sizeof (int) * 4 * nn * numStripes);
    for (int j=0 ; j<numStripes ; j++){
        lower[j] = (int*)((char *)(lower[0]) + j  * 4 * nn * sizeof(int)) ;
        upper[j] = (int*)((char *)(lower[j]) + 2 * nn * sizeof(int)) ;
    }
//    printf ("nn=%d %d\n", nn, sizeof (int) * 4 * nn * numStripes);
//    for (int j=0 ; j<numStripes ; j++)
//        printf("%d %ld %ld\n", j, (long long)lower[j] - (long long)lower[0], (long long) upper[j] - (long long)lower[0]);
    // printf("before parallel region!"); fflush(stdout);
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++){
        int imin = par * n / NUM_THREADS / 8 * 8;
        int imax = (par + 1) * n / NUM_THREADS / 8 * 8;
        if (par+1 == NUM_THREADS)
            imax = n;
        double h[d];
        for(int i=imin ; i<imax ; i++)
            lower[0][i+nn] = i/8 ;

        for(int j=1; j<numStripes ; j++){
            for(int a=0 ; a<d ; a++)
                h[a] = array[imin*d+a]+lowerDisplace[j][a] ;
            int a = imin ;
            int b = n-1 ;
            int m = (a+b)/2;
            while (b-a > 1){
                if(epsilonGridCompare(array + m*d, h) >= 0)
                    b = m ;
                else
                    a = m ;
                m = (a+b)/2 ;
            }
//            while (m < n && epsilonGridCompare(array + m * d, h) < 0)
//                m++;
//            printf("lower: imin = %d; m = %d %d %d\n", imin, a,m,b);
//            lower[j][imin+nn] = m/8 ;
            for(int i=imin ; i<imax ; i++) {
                for(int a=0 ; a<d ; a++)
                    h[a] = array[i*d+a]+lowerDisplace[j][a] ;
                while (m > 0 && epsilonGridCompare(array + m * d, h) >= 0)
                    m--;
                while (m < n && epsilonGridCompare(array + m * d, h) < 0)
                    m++;
                lower[j][i+nn] = m/8 ;
//                if(i == 17*1024){
//                    printf("lower j=%d, i=%d, m=%d\n", j, i, m);
//                    for (int a=0 ; a<d && a<3 ; a++)
//                        printf("%5.3lf ", array[m*d-d+a]);
//                    for (int a=0 ; a<d && a<d ; a++)
//                        printf("%3.1lf ", (double)(int)(array[m*d-d+a]/epsilon)*epsilon);
//                    printf("\n");
//                    for (int a=0 ; a<d && a<3 ; a++)
//                        printf("%5.3lf ", array[m*d+a]);
//                    for (int a=0 ; a<d && a<d ; a++)
//                        printf("%3.1lf ", (double)(int)(array[m*d+a]/epsilon)*epsilon);
//                    printf("\n");
//                }
            }
        }
// if (par == 0)
// {
//         printf("*");fflush(stdout);
// }
        for(int j=0; j<numStripes ; j++){
            for(int a=0 ; a<d ; a++)
                h[a] = array[imin*d+a]+upperDisplace[j][a] ;
            int a = imin ;
            int b = n-1 ;
            int m = (a+b)/2;
            while (b-a > 1){
                if(epsilonGridCompare(array + m*d, h) >= 0)
                    b = m ;
                else
                    a = m ;
                m = (a+b)/2 ;
            }
//            while (m < n && epsilonGridCompare(array + m * d, h) < 0)
//                m++;
//            printf("upper: imin = %d; m = %d %d %d\n", imin, a, m, b);
//            upper[j][imin+nn] = (m+7)/8 ;

    // printf("%d, ", omp_get_thread_num());fflush(stdout);

            for(int i=imin ; i<imax ; i++) {
                for(int a=0 ; a<d ; a++)
                    h[a] = array[i*d+a]+upperDisplace[j][a] ;
                while (m > 0 && epsilonGridCompare(array + m * d, h) >= 0)
                    m--;
                while (m < n && epsilonGridCompare(array + m * d, h) < 0)
                    m++;
                upper[j][i+nn] = (m+7)/8 ;
//                if(i == 17*1024){
//                    printf("upper j=%d, i=%d, m=%d\n", j, i, m);
//                    for (int a=0 ; a<d && a<3 ; a++)
//                        printf("%5.3lf ", array[m*d-d+a]);
//                    for (int a=0 ; a<d && a<d ; a++)
//                        printf("%3.1lf ", (double)(int)(array[m*d-d+a]/epsilon)*epsilon);
//                    printf("\n");
//                    for (int a=0 ; a<d && a<3 ; a++)
//                        printf("%5.3lf ", array[m*d+a]);
//                    for (int a=0 ; a<d && a<d ; a++)
//                        printf("%3.1lf ", (double)(int)(array[m*d+a]/epsilon)*epsilon);
//                    printf("\n");
//                }
            }
        }
        if(self) {
            double epsilon22 = epsilon * epsilon / 2;
            for (int i = imin; i < imax; i++) {
                double h = epsilon22;
                for (int j = 0; j < d; j++)
                    h -= array[i * d + j] * array[i * d + j];
                self[i] = h / 2;
            }
        }
        for(int j=0; j<numStripes ; j++){
            for (int i=imin/8 ; i<imax/8 ; i++)
                lower[j][i+nn/8] = min(min(min(lower[j][i*8+nn],lower[j][i*8+nn+1]),min(lower[j][i*8+nn+2], lower[j][i*8+nn+3])),
                        min(min(lower[j][i*8+nn+4],lower[j][i*8+nn+5]),min(lower[j][i*8+nn+6], lower[j][i*8+nn+7]))) ;
            for (int i=imin/8 ; i<imax/8 ; i++)
                upper[j][i+nn/8] = max(max(max(upper[j][i*8+nn],upper[j][i*8+nn+1]),max(upper[j][i*8+nn+2], upper[j][i*8+nn+3])),
                        max(max(upper[j][i*8+nn+4],upper[j][i*8+nn+5]),max(upper[j][i*8+nn+6], upper[j][i*8+nn+7]))) ;
        }
        for (int j=1 ; j<numStripes ; j++)
            for (int i=imin/8 ; i<imax/8 ; i++)
                upper[j-1][i+nn/8] = min(upper[j-1][i+nn/8], lower[j][i+nn/8]);
// if (par == 0)
// {
//         printf("-");fflush(stdout);
// }
    }
    // printf("\nend parallel region!\n"); fflush(stdout);
    for (int j=0 ; j<numStripes ; j++){
        epsilonGridCompleteListMin(nn/8, lower[j]);
        epsilonGridCompleteListMax(nn/8, upper[j]);
    }
//    printf("End prepare stripes\n");
//    for(int i=0; i<60; i++){
//        for(int j=0; j<numStripes; j++)
//            printf("%d %d ", lower[j][i+nn/65536], upper[j][i+nn/65536]);
//        printf("\n");
//    }
}

static inline long long _custom_mm512_reduce_add_epi64(__m512i a){
//    __m256i low  = _mm512_cvtepi64_epi32(a);
//    low = _mm256_hadd_epi32(low, low);
//    __m128i ulow = _mm_hadd_epi32(_mm256_castsi256_si128(low),_mm256_castsi256_si128(low));
//    return _mm_cvtsi128_si32(ulow) + _mm_extract_epi32(ulow, 1);
    __m256i b = _mm512_castsi512_si256(a) + _mm512_extracti64x4_epi64(a,1);
    __m128i c = _mm256_castsi256_si128(b) + _mm256_extracti128_si256(b,1);//_mm256_extracti64x2_epi64(b,1);
    return _mm_cvtsi128_si64(c) + _mm_extract_epi64(c, 1);
}

static inline void transpose_dx8(size_t n, size_t d, double *EGO_array) {
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++) {
        veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};
        veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};
        veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};
        veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};
        veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};
        veci64 mask6 = {12, 13, 1<4, 15, 4, 5, 6, 7};
        int ddd = (d+7)/8*8 ;
        double help[ddd*8]__attribute__((aligned(64)));
        for (int i = par * n / NUM_THREADS / 8; i < (par + 1) * n / NUM_THREADS / 8; i++) {
            for (int k = 0; k < d; k += 8) {
                vec v1 = _mm512_load_pd(EGO_array + 8 * (i*d + k));
                vec v2 = _mm512_load_pd(EGO_array + 8 * (i*d + k + 1));
                vec v3 = _mm512_load_pd(EGO_array + 8 * (i*d + k + 2));
                vec v4 = _mm512_load_pd(EGO_array + 8 * (i*d + k + 3));
                vec v5 = _mm512_load_pd(EGO_array + 8 * (i*d + k + 4));
                vec v6 = _mm512_load_pd(EGO_array + 8 * (i*d + k + 5));
                vec v7 = _mm512_load_pd(EGO_array + 8 * (i*d + k + 6));
                vec v8 = _mm512_load_pd(EGO_array + 8 * (i*d + k + 7));
                transposeAVX512(v1, v2, v3, v4, v5, v6, v7, v8);
                _mm512_store_pd(help + k, v1);
                _mm512_store_pd(help + k+ddd, v2);
                _mm512_store_pd(help + k+2*ddd, v3);
                _mm512_store_pd(help + k+3*ddd, v4);
                _mm512_store_pd(help + k+4*ddd, v5);
                _mm512_store_pd(help + k+5*ddd, v6);
                _mm512_store_pd(help + k+6*ddd, v7);
                _mm512_store_pd(help + k+7*ddd, v8);
            }
            for (int k=0 ; k<8 ; k++)
                memcpy(EGO_array + (i*8+k)*d, help + ddd*k, d * sizeof (double));
        }
    }
}


int cmp_reorder_dim(const void * a, const void *b){
    long long A = costref[*(int*) a];
    long long B = costref[*(int*) b];
    if(A<B) return -1;
    if(A>B) return 1;
    return 0;
}


void outputStatistics(int n, int d, double epsilon, double *array, int *reorder_dim){
    // printf("Statistics\n");
    //    int first_int[d], last_int[d];
    //    for (int j=0 ; j<d ; j++)
    //        first_int[j] = last_int[j] = (int) floor(array[j]/epsilon);
    //    for (int i=0 ; i<n ; i++)
    //        for (int j=0 ; j<d ; j++){
    //            first_int[j] = min(first_int[j], (int) floor(array[i*d+j]/epsilon));
    //            last_int[j] = max(last_int[j], (int) floor(array[i*d+j]/epsilon));
    //        }
    costref = (long long *) malloc(d*sizeof(long long));
    for (int j = 0; j < d; j++) {
        double first_int = array[j];
        double last_int = first_int;
        int testcount = 0;
#pragma omp parallel for reduction(min:first_int) reduction(max:last_int) reduction(+:testcount)
        for(int par=0 ; par<NUM_THREADS ; par++)
        for (int i = par+1; i < n; i+=(int)(-log(1.-drand48())*100.*NUM_THREADS)){
            first_int = min(first_int, array[i * d + j]);
            last_int = max(last_int, array[i * d + j]);
            testcount ++;
        }
        first_int = floor(first_int/epsilon);
        last_int = floor(last_int/epsilon);
        if (first_int < -100000 || last_int > 100000) {
#ifdef OUTPUT
            printf("%3d %e %e ------------------ (VALUE RANGE EXCEEDED)\n", j, first_int, last_int);
#endif
        } else {
            int size = (int) last_int - (int) first_int + 1;
            long long * histo = (long long *) calloc(size, sizeof (long long));
            //memset(histo, 0, size*sizeof(long long));
            int testcount2 = 0;
#pragma omp parallel for reduction(+:histo[:size]) reduction(+:testcount2)
            for (int par = 0; par < NUM_THREADS; par++)
                for (int i = par + 1; i < n; i += (int) (-log(1. - drand48())*100. * NUM_THREADS), testcount2++)
                    histo[max(0, min(size - 1, (int) floor(array[i * d + j] / epsilon) - (int) first_int))]++;
            long long ref = 0;
            for (int i = 0; i < size; i++)
                ref += histo[i]*(histo[i] - 1) / 2;
            for (int i = 1; i < size; i++)
                ref += histo[i - 1] * histo[i];
            costref[j] = ref;
#ifdef OUTPUT
            printf("%3d %8d %8d %8d %20lld   [%lld", j, (int) first_int, (int) last_int, size, ref, histo[0]);

            for (int i = 1; i < size && i < 10; i++)
                printf(", %lld", histo[i]);
            if (size > 10)
                printf(", ...] %d %d\n", testcount, testcount2);
            else
                printf("]\n");
#endif
            free(histo);
        }
    }
    int * reorder_rev = (int*) malloc ((d+8)*sizeof(int));
    for (int j=0 ; j<d+8 ; j++)
        reorder_rev[j] = j;
    qsort(reorder_rev, d, sizeof(int), cmp_reorder_dim);
#ifdef OUTPUT
    for (int j=0 ; j<d+8 ; j++)
        printf("%2d %2d %lld\n", j, reorder_rev[j], j<d?costref[reorder_rev[j]]:0);
#endif
    // reorder_dim = (int*) malloc ((d+8)*sizeof(int));
    for (int j=0 ; j<d+8 ; j++)
        reorder_dim[reorder_rev[j]] = j;
    free(reorder_rev);
    free(costref);
}

void omp_qsort(void* l, size_t num, size_t size, int (*compar)(const void*,const void*)){
//    printf ("Qsort %ld %d %d %d\n", (long long)l, num, omp_get_num_threads(), omp_get_thread_num());
    if(num <= min_size_qsort){
        qsort(l, num, size, compar);
        return ;
    }
    // a < b < c || c < b < a --> b
    // b < a < c || c < a < b --> a
    // a < c < b || b < c < a --> c
    char b[size] ;
    char *v = (char*)l + (num-1)*size ;
    char *r = v ;
    char *i = (char*)l - size ;
    char *j = (char*)l+num/2*size ;
    if(compar(l, j) < 0){
        if(compar(j, r) < 0){
            hilbert_swap(j,r)
        } else {
            if(compar(r, l) < 0)
                hilbert_swap(l,r)
            }
    } else {
        if(compar(r,j) < 0){
            hilbert_swap(j,r)
        } else if(compar(l, r) < 0)
            hilbert_swap(l,r)
    }
    j = r ;
    for (;;) {
        while (compar(i+=size, v) < 0) ;
        while (compar(v, j-=size) < 0)
            if(j==i) break ;
        if(i>=j) break ;
        hilbert_swap(i,j)
    }
    hilbert_swap(i,r)
#pragma omp task
        omp_qsort (l, ((long long)i-(long long)l) / size, size, compar);
#pragma omp task
        omp_qsort (i+size, num - ((long long)i-(long long)l)/size - 1, size, compar);
}

void reorder_dimensions(int n, int d, double *array, int *reorder_dim) {
    if (d > 128) {
#pragma omp parallel for
        for (int par = 0; par < NUM_THREADS; par++) {
            double help[d + 8]__attribute__((aligned(64)));
            for (int i = par * n / NUM_THREADS; i < (par + 1) * n / NUM_THREADS; i++) {
                memcpy(help, array + i*d, d * sizeof (double));
                for (int j = 0; j < d; j++)
                    array[reorder_dim[j]+i*d] = help[j];
            }
        }
        return;
    }
    // printf("first\n");fflush(stdout);
    // printf("reorder_dim[0]:%d\n", reorder_dim[0]);fflush(stdout);
    // printf("second");fflush(stdout);
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++) {
        __m256i i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15;
        __m256i * rd = (__m256i *) reorder_dim;
        switch ((d - 1) / 8) {
            case 15: i15 = _mm256_loadu_si256(rd++); //rd += 8;
            case 14: i14 = _mm256_loadu_si256(rd++); //rd += 8;
            case 13: i13 = _mm256_loadu_si256(rd++); //rd += 8;
            case 12: i12 = _mm256_loadu_si256(rd++); //rd += 8;
            case 11: i11 = _mm256_loadu_si256(rd++); //rd += 8;
            case 10: i10 = _mm256_loadu_si256(rd++); //rd += 8;
            case 9: i9 = _mm256_loadu_si256(rd++); //rd += 8;
            case 8: i8 = _mm256_loadu_si256(rd++); //rd += 8;
            case 7: i7 = _mm256_loadu_si256(rd++); //rd += 8;
            case 6: i6 = _mm256_loadu_si256(rd++); //rd += 8;
            case 5: i5 = _mm256_loadu_si256(rd++); //rd += 8;
            case 4: i4 = _mm256_loadu_si256(rd++); //rd += 8;
            case 3: i3 = _mm256_loadu_si256(rd++); //rd += 8;
            case 2: i2 = _mm256_loadu_si256(rd++); //rd += 8;
            case 1: i1 = _mm256_loadu_si256(rd++); //rd += 8;
            case 0: i0 = _mm256_loadu_si256(rd);
        }
        __mmask8 k = 256 - (1 << ((128 - d) % 8));
        vec x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15;
        for (int i = par * n / NUM_THREADS; i < (par + 1) * n / NUM_THREADS; i++) {
            double * xd = array + i*d;
            switch ((d - 1) / 8) {
                case 15: x15 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 14: x14 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 13: x13 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 12: x12 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 11: x11 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 10: x10 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 9: x9 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 8: x8 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 7: x7 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 6: x6 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 5: x5 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 4: x4 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 3: x3 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 2: x2 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 1: x1 = _mm512_loadu_pd(xd);
                    xd += 8;
                case 0: x0 = _mm512_loadu_pd(xd);
            }
            xd = array + i*d;
            switch ((d - 1) / 8) {
                case 15: _mm512_i32scatter_pd(xd, i15, x15, 8);
                case 14: _mm512_i32scatter_pd(xd, i14, x14, 8);
                case 13: _mm512_i32scatter_pd(xd, i13, x13, 8);
                case 12: _mm512_i32scatter_pd(xd, i12, x12, 8);
                case 11: _mm512_i32scatter_pd(xd, i11, x11, 8);
                case 10: _mm512_i32scatter_pd(xd, i10, x10, 8);
                case 9: _mm512_i32scatter_pd(xd, i9, x9, 8);
                case 8: _mm512_i32scatter_pd(xd, i8, x8, 8);
                case 7: _mm512_i32scatter_pd(xd, i7, x7, 8);
                case 6: _mm512_i32scatter_pd(xd, i6, x6, 8);
                case 5: _mm512_i32scatter_pd(xd, i5, x5, 8);
                case 4: _mm512_i32scatter_pd(xd, i4, x4, 8);
                case 3: _mm512_i32scatter_pd(xd, i3, x3, 8);
                case 2: _mm512_i32scatter_pd(xd, i2, x2, 8);
                case 1: _mm512_i32scatter_pd(xd, i1, x1, 8);
                case 0: _mm512_mask_i32scatter_pd(xd, k, i0, x0, 8);
            }
        }
    }
}

// void test_ego_loop3_macro(size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, size_t *countresult, int stripes, int KBLOCK, double *sorttime){
//     size_t result = 0;
//     CUtilTimer sortTimer;
//     sortTimer.start();
//     EGO_PARALLEL_TRAN(n, d, epsilon, 5, array,NUM_THREADS)
//         sortTimer.stop();
//         // printf("Ego_n= %d\n", n);
//         // printf("ceilpowtwo(n)= %d\n", ceilpowtwo(EGO_n));
//         // for ( int __i = 0 ; __i < 2 ; __i++ ){
//         //     for ( int __j =0 ; __j < 5 ; __j++ ){
//         //         printf("lower[%d][%d]:%7d ", __i, __j, lower[__i][__j]);
//         //     }
//         //     printf("\n");
//         //     for ( int __j =0 ; __j < 5 ; __j++ ){
//         //         printf("upper[%d][%d]:%7d ", __i, __j, upper[__i][__j]);
//         //     }
//         //     printf("\n\n");
//         // }
//         #pragma omp parallel for reduction(+:result)
//     EGO_PREPARE
//         veci64 resultvec = _mm512_setzero_si512();
//         veci64 eights = _mm512_set1_epi64(8ll) ;
//     EGO_LOOP_TRAN{
//         resultvec += eights - _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63);
//     }
//     EGO_CONSOLIDATE{
//         result += _mm512_reduce_add_epi64(resultvec);
// //        double testres[8] __attribute__((aligned(64)));
// //        _mm512_store_epi64(testres, resultvec);
// //        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
//     }
//     EGO_END_TRAN
//     *sorttime = sortTimer.get_time();
//     *countresult = result;
// }

void test_ego_loop3_macro(size_t n, size_t d, double epsilon, double *array, size_t *countresult, int stripes, double *sortTime, double *indextime, double *loadpercent){
    CUtilTimer index_timer;
    long long result = 0;
    long long refinements = 0;
    unsigned long long savedload[5*NUM_THREADS];
    index_timer.start();
    EGO_PARALLEL_TRAN(n, d, epsilon, stripes, array)
        index_timer.stop();
        // printf("saving result...\n");
        // char filename[50];
        // sprintf (filename, "uniform_1000000x32_normalized_sorted_%f.bin", epsilon);
        // save_binary_file(array, n, d, filename);
        // printf("saved.\n");
        CUtilTimer total_timer; total_timer.start();
        *indextime = index_timer.get_time() - sortTimer.get_time();
        *sortTime = sortTimer.get_time();
        // printf("timestamp index ready %6.2f\n",timestamp()-starttimestamp);
#ifdef OUTPUT
        printf("overall_load: %lld / %lld (=n*(n-1)/2 / 64) ==> %f\n", overall_load, (long long)n/128*(n-1), (double)overall_load/n/(n-1)*128);
#endif
        *loadpercent = (double)overall_load/n/(n-1)*128;
        #pragma omp parallel for reduction(+:result) reduction(+:refinements)
    EGO_PREPARE
        veci64 resultvec = _mm512_setzero_si512();
        veci64 eights = _mm512_set1_epi64(8ll) ;
        long long refineload = 0;
    EGO_LOOP_TRAN {
        resultvec += eights - _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63) ;
        refineload ++;
    }
    EGO_CONSOLIDATE{
        result += _custom_mm512_reduce_add_epi64(resultvec);
        refinements += refineload ;
        int curload=0;
        for(int i=loadstart[par] ; i<loadstart[par+1] ; i++)
            for(int s=0 ; s<EGO_stripes ; s++)
                curload += upper[s][i+nn/8] - lower[s][i+nn/8];
        total_timer.stop();
#ifdef OUTPUT
        printf("Consolidate %6.2f %d %d %d %d %d %lld %lld\n",total_timer.get_time(), par, omp_get_thread_num(), loadstart[par], loadstart[par+1]-loadstart[par], curload, refineload, result);
#endif

//        double testres[8] __attribute__((aligned(64)));
//        _mm512_store_epi64(testres, resultvec);
//        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
    }
    EGO_END_TRAN
    // printf("result %ld\n", result);
    *countresult = result;

    // printf("refinements %ld Mio (%ld 8x8)\n", refinements*64/1000000, refinements);
//    for(int par=0 ; par<NUM_THREADS ; par++, printf("\n"))
//        for(int s=0 ; s<5 ; s++)
//            printf("%ld ",savedload[NUM_THREADS*s+par]);
}

void test_ego_loop3_noself(const size_t nA, const size_t nB, const int d, const double epsilon, double *arrayA, double *arrayB, size_t *countresult, const int activedims, double *sortTime, double *indextime, double *loadpercent){
    CUtilTimer index_timer,total_timer,sortTimer;
    long long result = 0;
    long long refinements = 0;
    unsigned long long savedload[5*NUM_THREADS];
    double starttimestamp = timestamp() ;
    index_timer.start();
    EGO_PARALLEL_TRAN_NOSELF(nA, nB, d, epsilon, 2, arrayA, arrayB)
        index_timer.stop();
        // printf("timestamp index ready %6.2f\n",timestamp()-starttimestamp);
//        for(int i=0 ; i<NUM_THREADS + 4 ; i+=4)
//            printf("%2d %9d %9d %9d %9d\n", i, loadstart[i], loadstart[i+1], loadstart[i+2], loadstart[i+3]);
//        for (int x=0 ; x<NUM_THREADS ; x++){
//            long long cum_load = 0;
//            for(int i = loadstart[x] ; i<loadstart[x+1] ; i++)
//                for (int j = 0; j < EGO_stripes; j++)
//                    cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];
//            if (x%4 == 0) printf("\n %2d", x);
//            printf("%9ld (%5.2f%%)  ", cum_load, 100.*(double)cum_load/(double)overall_load);
//        }
        // printf("overall_load: %ld / %ld (=n/8*(n/8+1)/2) ==> %9.6f %%\n", overall_load, (long long)nA/8*(nA/8+1), (double)overall_load/(nA/8)/(nA/8-1)*200);
       total_timer.start();
       *indextime = index_timer.get_time() - sortTimer.get_time();
       *sortTime = sortTimer.get_time();

       *loadpercent = ((double)overall_load/(nA/8)/(nA/8-1)*200);
       #pragma omp parallel for proc_bind(close) reduction(+:result) reduction(+:refinements)
    EGO_PREPARE
        veci64 resultvec = _mm512_setzero_si512();
        veci64 eights = _mm512_set1_epi64(8ll) ;
        long long refineload = 0;
    EGO_LOOP_TRAN_NOSELF {
        resultvec += eights - _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63)
                            - _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63) ;
        refineload ++;
    }
    EGO_CONSOLIDATE{
        result += _custom_mm512_reduce_add_epi64(resultvec);
        refinements += refineload ;
        int curload=0;
        for(int i=loadstart[par] ; i<loadstart[par+1] ; i++)
            for(int s=0 ; s<EGO_stripes ; s++)
                curload += upper[s][i+nn/8] - lower[s][i+nn/8];
        total_timer.stop();
#ifdef OUTPUT
        printf("Consolidate %6.2f %d %d %d %d %d %lld %lld\n",timestamp()-starttimestamp, par, omp_get_thread_num(), loadstart[par], loadstart[par+1]-loadstart[par], curload, refineload, result);
#endif

//        double testres[8] __attribute__((aligned(64)));
//        _mm512_store_epi64(testres, resultvec);
//        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
    }
    EGO_END_TRAN_NOSELF

    *countresult = result;
    // printf("result %ld\n", result);
    // printf("refinements %ld Mio (%ld 8x8)\n", refinements*64/1000000, refinements);

//    for(int par=0 ; par<NUM_THREADS ; par++, printf("\n"))
//        for(int s=0 ; s<5 ; s++)
//            printf("%ld ",savedload[NUM_THREADS*s+par]);
}


// void test_ego_loop3_macro_queue(size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, size_t *countresult, int stripes, int KBLOCK, boost::lockfree::queue<join_pair> &jpartners, double *sorttime){
//     size_t result = 0;
//     CUtilTimer sortTimer;
//     sortTimer.start();
//     EGO_PARALLEL_TRAN(n, d, epsilon, 5, array,NUM_THREADS)
//         sortTimer.stop();
//         #pragma omp parallel
//         {
//             boost::lockfree::queue<join_pair> private_queue(128);
//             #pragma omp for
//             EGO_PREPARE
//                 veci64 resultvec = _mm512_setzero_si512();
//                 veci64 eights = _mm512_set1_epi64(8ll) ;
//             EGO_LOOP_TRAN{
//                 resultvec += eights - _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63)
//                                     - _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63)
//                                     - _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63)
//                                     - _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63)
//                                     - _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63)
//                                     - _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63)
//                                     - _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63)
//                                     - _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63);
//                 double outarray[64] __attribute__((aligned(64)));
//                 _mm512_store_pd(&outarray[0], sum1);
//                 _mm512_store_pd(&outarray[8], sum2);
//                 _mm512_store_pd(&outarray[16], sum3);
//                 _mm512_store_pd(&outarray[24], sum4);
//                 _mm512_store_pd(&outarray[32], sum5);
//                 _mm512_store_pd(&outarray[40], sum6);
//                 _mm512_store_pd(&outarray[48], sum7);
//                 _mm512_store_pd(&outarray[54], sum8);
//
//                 for ( int ii = 0 ; ii < 64 ; i++ ){
//                     if ( outarray[ii] > 0.0 ){
//                         join_pair p;
//                         p.p1 = i*8 + (ii % 8);  // i index
//                         p.p2 = j*8 + (ii / 8); // j index
//                         while ( ! private_queue.push(p) );
//                     }
//                 }
//             }
//             EGO_CONSOLIDATE{
//                 result += _mm512_reduce_add_epi64(resultvec);
//         //        double testres[8] __attribute__((aligned(64)));
//         //        _mm512_store_epi64(testres, resultvec);
//         //        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
//             } // END_CONSOLIDATE
//         }// END omp parallel
//     EGO_END_TRAN
//     *sorttime = sortTimer.get_time();
//     *countresult = result;
// }

// void test_ego_loop3_long_queue(size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, size_t *countresult, int stripes, int KBLOCK, boost::lockfree::queue<join_pair> &jpartners){
//     size_t result = 0;
//
//     size_t EGO_n = (n);
//     EGO_d = (d);
//     EGO_epsilon = (epsilon);
//     double * EGO_array = (array);
//     int EGO_blocks = (EGO_d + KBLOCK - 1) / KBLOCK;
//     int EGO_stripes = (stripes);
//     epsilonGridOrdering(EGO_n, EGO_d, EGO_epsilon, EGO_array,NUM_THREADS);
//     int nn = ceilpowtwo(EGO_n);
//     int **lower = (int **) ddr_alloc (EGO_stripes*sizeof(int*));
//     int **upper = (int **) ddr_alloc (EGO_stripes*sizeof(int*));
//     double *self = callocA64(sizeof (double) * EGO_n * EGO_blocks);
//     prepareStripes(EGO_n, EGO_d, EGO_stripes, EGO_epsilon, EGO_array, lower, upper, (double *)0, NUM_THREADS);
//     EGO_epsilon = EGO_epsilon * EGO_epsilon / 2;
//     for (int i=0 ; i<EGO_n ; i++){
//         double h=EGO_epsilon;
//         for(int j=0 ; j<2*KBLOCK && j<EGO_d ; j++)
//             h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];
//         self[i/8*8*EGO_blocks+i%8]=h/2;
//     }
//     for(int k=2 ; k<EGO_blocks ; k++)
//         for (int i=0 ; i<EGO_n ; i++){
//             double h=0;
//             for(int j=k * KBLOCK ; j<(k+1) * KBLOCK && j<EGO_d ; j++)
//             h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];
//         self[(i/8*EGO_blocks+k)*8+i%8]=h/2;
//     }
//     transpose_8xd(EGO_n, EGO_d, EGO_array,NUM_THREADS);
//     long long overall_load = 0;
//     for (int i = 0; i < EGO_n / 8; i++)
//         for (int j = 0; j < EGO_stripes; j++)
//             overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];
//     int loadstart[NUM_THREADS + 2];
//     for (int i = 0; i <= NUM_THREADS; i++)
//         loadstart[i] = 0;
//     long long cum_load = 0;
//     for (int i = 0; i < EGO_n / 8; i++) {
//         for (int j = 0; j < EGO_stripes; j++)
//             cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];
//         loadstart[cum_load * NUM_THREADS / overall_load + 1] = i;
//     }
//     loadstart[NUM_THREADS] = EGO_n / 8;
//     for (int i = 1; i <= NUM_THREADS; i++)
//         if (loadstart[i] == 0)
//             loadstart[i] = loadstart[i - 1];
//
//         #pragma omp parallel for reduction(+:result)
//         for (int par = 0; par < NUM_THREADS; par++) {
//
//         veci64 resultvec = _mm512_setzero_si512();
//         veci64 eights = _mm512_set1_epi64(8ll) ;
//
//         int imin = loadstart[par];
//         int imax = loadstart[par + 1];
//         for (int s = 0; s < EGO_stripes; s++) {
//             int i = 0;
//             int j = 0;
//             register veci64 const1 = _mm512_set1_epi64(1LL) ;
//             register veci64 const2 = _mm512_add_epi64(const1, const1) ;
//             register veci64 const3 = _mm512_add_epi64(const1, const2) ;
//             register veci64 const4 = _mm512_add_epi64(const1, const3) ;
//             register veci64 const5 = _mm512_add_epi64(const1, const4) ;
//             register veci64 const6 = _mm512_add_epi64(const1, const5) ;
//             register veci64 const7 = _mm512_add_epi64(const1, const6) ;
//             register veci64 const0 = _mm512_setzero_si512();
//             FGF_HILBERT_FOR(i, j, EGO_n / 8, EGO_n / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],
//                     FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&
//                     FURHIL_ub1 - 1 >= lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&
//                     FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {
//                 register vec vi = _mm512_load_pd(self + i * 8 * EGO_blocks);
//                 register vec vj = _mm512_load_pd(self + j * 8 * EGO_blocks);
//                 register vec sum1 = vi + _mm512_permutexvar_pd(const0, vj);
//                 register vec sum2 = vi + _mm512_permutexvar_pd(const1, vj);
//                 register vec sum3 = vi + _mm512_permutexvar_pd(const2, vj);
//                 register vec sum4 = vi + _mm512_permutexvar_pd(const3, vj);
//                 register vec sum5 = vi + _mm512_permutexvar_pd(const4, vj);
//                 register vec sum6 = vi + _mm512_permutexvar_pd(const5, vj);
//                 register vec sum7 = vi + _mm512_permutexvar_pd(const6, vj);
//                 register vec sum8 = vi + _mm512_permutexvar_pd(const7, vj);
//                 vi = _mm512_load_pd(EGO_array + (EGO_d * i) * 8);
//                 vj = _mm512_load_pd(EGO_array + (EGO_d * j) * 8);
//                 sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);
//                 sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);
//                 sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);
//                 sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);
//                 sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);
//                 sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);
//                 sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);
//                 sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);
//                 int k; for (k=1 ; k<d ; k++){
//                     if(k % KBLOCK == 0 && k > KBLOCK){
//                         register veci64 allind = _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63);
//                         allind += _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63);
//                         allind += _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63);
//                         allind += _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63);
//                         allind += _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63);
//                         allind += _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63);
//                         allind += _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63);
//                         allind += _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63);
//                         if(_mm512_reduce_add_epi64(allind) >= 64) {k=d+1; break;}
//                         vi = _mm512_load_pd(self + (i * EGO_blocks + k/KBLOCK) * 8);
//                         vj = _mm512_load_pd(self + (j * EGO_blocks + k/KBLOCK) * 8);
//                         sum1 += vi + _mm512_permutexvar_pd(const0, vj);
//                         sum2 += vi + _mm512_permutexvar_pd(const1, vj);
//                         sum3 += vi + _mm512_permutexvar_pd(const2, vj);
//                         sum4 += vi + _mm512_permutexvar_pd(const3, vj);
//                         sum5 += vi + _mm512_permutexvar_pd(const4, vj);
//                         sum6 += vi + _mm512_permutexvar_pd(const5, vj);
//                         sum7 += vi + _mm512_permutexvar_pd(const6, vj);
//                         sum8 += vi + _mm512_permutexvar_pd(const7, vj);
//                     }
//                     vi = _mm512_load_pd(EGO_array + (EGO_d * i + k) * 8);
//                     vj = _mm512_load_pd(EGO_array + (EGO_d * j + k) * 8);
//                     sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);
//                     sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);
//                     sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);
//                     sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);
//                     sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);
//                     sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);
//                     sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);
//                     sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);
//                 }
//                 if(k<=d){{
//                     if (i==j){
//                         sum1 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum1), 255, 0xbff0000000000000ull));
//                         sum2 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum2), 254, 0xbff0000000000000ull));
//                         sum3 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum3), 252, 0xbff0000000000000ull));
//                         sum4 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum4), 248, 0xbff0000000000000ull));
//                         sum5 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum5), 240, 0xbff0000000000000ull));
//                         sum6 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum6), 224, 0xbff0000000000000ull));
//                         sum7 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum7), 192, 0xbff0000000000000ull));
//                         sum8 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum8), 128, 0xbff0000000000000ull));
//                     }
//
//     {
//         resultvec += eights - _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63)
//                             - _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63) ;
//     }
//
//
//     }}} FGF_HILBERT_END(i, j); }
//
//         {
//             result += _mm512_reduce_add_epi64(resultvec);
//     //        double testres[8] __attribute__((aligned(64)));
//     //        _mm512_store_epi64(testres, resultvec);
//     //        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
//         }
//
//     }transpose_dx8(EGO_n, EGO_d, EGO_array, NUM_THREADS);}
//
//     // *countresult = result;
// }


void prepareStripesNoSelf(int nA, int nB, int d, int activeDimensions, double epsilon, double *A, double *B, int ** lower, int **upper, double *selfA, double *selfB) {
    double starttimestamp = timestamp();
    int numStripes = 1;
    for (int i = 0; i < activeDimensions; i++)
        numStripes *= 3;
    double lowerDisplace[numStripes][d];
    double upperDisplace[numStripes][d];
    for (int i = 0; i < numStripes; i++)
        for (int j = 0; j < d; j++) {
            lowerDisplace[i][j] = (-epsilon);
            upperDisplace[i][j] = epsilon;
        }
    for (int i = 0; i < numStripes; i++) {
        int power = numStripes / 3;
        for (int j = 0; j < activeDimensions; j++) {
            lowerDisplace[i][j] = upperDisplace[i][j] = epsilon * (double) (i / power % 3 - 1);
            power /= 3;
        }
    }
    int nn = ceilpowtwo(nA);
    lower[0] = (int *) callocA64(sizeof (int) * 4 * nn * numStripes);
    for (int j=0 ; j<numStripes ; j++){
        lower[j] = (int*)((char *)(lower[0]) + j  * 4 * nn * sizeof(int)) ;
        upper[j] = (int*)((char *)(lower[j]) + 2 * nn * sizeof(int)) ;
    }
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++){
        int imin = par * nA / NUM_THREADS / 8 * 8;
        int imax = (par + 1) * nA / NUM_THREADS / 8 * 8;
        if (par+1 == NUM_THREADS)
            imax = nA;
        double h[d];
        for(int j=0; j<numStripes ; j++){
            for(int a=0 ; a<d ; a++)
                h[a] = A[imin*d+a]+lowerDisplace[j][a] ;
            int a = 0 ;
            int b = nB-1 ;
            int m = (a+b)/2;
            while (b-a > 1){
                if(epsilonGridCompare(B + m*d, h) >= 0)
                    b = m ;
                else
                    a = m ;
                m = (a+b)/2 ;
            }
            for(int i=imin ; i<imax ; i++) {
                for(int a=0 ; a<d ; a++)
                    h[a] = A[i*d+a]+lowerDisplace[j][a] ;
                while (m > 0 && epsilonGridCompare(B + m * d, h) >= 0)
                    m--;
                while (m < nB && epsilonGridCompare(B + m * d, h) < 0)
                    m++;
                lower[j][i+nn] = m/8 ;
            }
        }
        for(int j=0; j<numStripes ; j++){
            for(int a=0 ; a<d ; a++)
                h[a] = A[imin*d+a]+upperDisplace[j][a] ;
            int a = imin ;
            int b = nB-1 ;
            int m = (a+b)/2;
            while (b-a > 1){
                if(epsilonGridCompare(B + m*d, h) >= 0)
                    b = m ;
                else
                    a = m ;
                m = (a+b)/2 ;
            }
            for(int i=imin ; i<imax ; i++) {
                for(int a=0 ; a<d ; a++)
                    h[a] = A[i*d+a]+upperDisplace[j][a] ;
                while (m > 0 && epsilonGridCompare(B + m * d, h) >= 0)
                    m--;
                while (m < nB && epsilonGridCompare(B + m * d, h) < 0)
                    m++;
                upper[j][i+nn] = (m+7)/8 ;
            }
        }
        double epsilon22 = epsilon * epsilon / 2;
        if(selfA) {
            for (int i = imin; i < imax; i++) {
                double h = epsilon22;
                for (int j = 0; j < d; j++)
                    h -= A[i * d + j] * A[i * d + j];
                selfA[i] = h / 2;
            }
        }
        if(selfB) {
            for (int i = par*nB/NUM_THREADS ; i < (par+1)*nB/NUM_THREADS; i++) {
                double h = epsilon22;
                for (int j = 0; j < d; j++)
                    h -= B[i * d + j] * B[i * d + j];
                selfB[i] = h / 2;
            }
        }
        for(int j=0; j<numStripes ; j++){
            for (int i=imin/8 ; i<imax/8 ; i++)
                lower[j][i+nn/8] = min(min(min(lower[j][i*8+nn],lower[j][i*8+nn+1]),min(lower[j][i*8+nn+2], lower[j][i*8+nn+3])),
                        min(min(lower[j][i*8+nn+4],lower[j][i*8+nn+5]),min(lower[j][i*8+nn+6], lower[j][i*8+nn+7]))) ;
            for (int i=imin/8 ; i<imax/8 ; i++)
                upper[j][i+nn/8] = max(max(max(upper[j][i*8+nn],upper[j][i*8+nn+1]),max(upper[j][i*8+nn+2], upper[j][i*8+nn+3])),
                        max(max(upper[j][i*8+nn+4],upper[j][i*8+nn+5]),max(upper[j][i*8+nn+6], upper[j][i*8+nn+7]))) ;
        }
        for (int j=1 ; j<numStripes ; j++)
            for (int i=imin/8 ; i<imax/8 ; i++)
                upper[j-1][i+nn/8] = min(upper[j-1][i+nn/8], lower[j][i+nn/8]);
#ifdef OUTPUT
        printf("Thread %d ready %7.2f %d %d\n", par, timestamp()-starttimestamp, imin, imax);
#endif
    }
    for (int j=0 ; j<numStripes ; j++){
        epsilonGridCompleteListMin(nn/8, lower[j]);
        epsilonGridCompleteListMax(nn/8, upper[j]);
    }
}
