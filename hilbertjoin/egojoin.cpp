
#include "egojoin.h"

double EGO_epsilon;
int EGO_d;
int min_size_qsort;

int epsilonGridCompare(const void *a, const void *b) {
    double * A = (double *) a;
    double * B = (double *) b;
    for (int i = 0; i < EGO_d; i++) {
        int h = (int) (floor(A[i] / EGO_epsilon) - floor(B[i] / EGO_epsilon));
        if (h != 0) return h;
    }
    return 0;
}


int test_ego_loop(size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, long long *result){
    int KBLOCK=8;
    *result = 0l;
    long long iresult = 0l;
    EGO_PARALLEL(n, d, epsilon, array, NUM_THREADS)
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
    printf("count of join partners: %d\n", *iresult);
#endif
    *result = iresult;

}

int test_ego_loop3(size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, long long *result){
    *result = 0l;
    int KBLOCK=8;
    long long iresult=0l;
    omp_set_num_threads(NUM_THREADS);
    EGO_PARALLEL_TRAN(n, d, epsilon, 5, array, NUM_THREADS)
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
        iresult += _mm512_reduce_add_epi64(resultvec);
//        double testres[8] __attribute__((aligned(64)));
//        _mm512_store_epi64(testres, resultvec);
//        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
    }
    EGO_END_TRAN
#ifdef OUTPUT
    printf("result %d\n", iresult);
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

void epsilonGridOrdering(size_t n, size_t d, double epsilon, double* array, size_t NUM_THREADS) {
    EGO_epsilon = epsilon;
    EGO_d = d;
    min_size_qsort = n/NUM_THREADS*2 ;
#pragma omp parallel
#pragma omp master
    omp_qsort(array, n, d * sizeof (double), epsilonGridCompare);
#pragma omp taskwait
//    for (int i=1 ; i<n ; i++)
//        if(epsilonGridCompare(array+d*(i-1), array+d*i)>0)
//            printf("Qsort Error %ld\n", (long long) (array+d*i));
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
        list[i] = max(list[2 * i], list[2 * i + 1]);
}

void epsilonGridCompleteListMin(size_t n, int *list) {
    int m = ceilpowtwo(n);
    for (size_t i = n; i < m; i++)
        list[i+m] = 0;
    if (n % 2)
        list[n] = list[2 * n - 1];
    for (size_t i = n - 1; i > 0; i--)
        list[i] = min(list[2 * i], list[2 * i + 1]);
}

static inline void transpose_8xd(size_t n, size_t d, double *EGO_array, size_t NUM_THREADS) {
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


void prepareStripes(size_t n, size_t d, int numStripes, double epsilon, double *array, int ** lower, int ** upper, double *self,size_t NUM_THREADS){
//                    for (int a=0 ; a<d && a<20 ; a++)
//                        printf("%5.3lf ", array[17*1024*d+a]);
//                    printf("\n");
//                    for (int a=0 ; a<d && a<d ; a++)
//                        printf("%3.1lf ", (double)(int)(array[17*1024*d+a]/epsilon)*epsilon);
//                    printf("\n");
    double lowerDisplace[numStripes][d];
    double upperDisplace[numStripes][d];
    EGO_epsilon = epsilon;
    EGO_d = d;
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
    }
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

static inline long long _mm512_reduce_add_epi64(__m512i a){
//    __m256i low  = _mm512_cvtepi64_epi32(a);
//    low = _mm256_hadd_epi32(low, low);
//    __m128i ulow = _mm_hadd_epi32(_mm256_castsi256_si128(low),_mm256_castsi256_si128(low));
//    return _mm_cvtsi128_si32(ulow) + _mm_extract_epi32(ulow, 1);
    __m256i b = _mm512_castsi512_si256(a) + _mm512_extracti64x4_epi64(a,1);
    __m128i c = _mm256_castsi256_si128(b) + _mm256_extracti128_si256(b,1);//_mm256_extracti64x2_epi64(b,1);
    return _mm_cvtsi128_si64(c) + _mm_extract_epi64(c, 1);
}

static inline void transpose_dx8(size_t n, size_t d, double *EGO_array, size_t NUM_THREADS) {
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++) {
        veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};
        veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};
        veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};
        veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};
        veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};
        veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};
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

void omp_qsort (void* l, size_t num, size_t size, int (*compar)(const void*,const void*)){
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
            swap(j,r);
        } else {
            if(compar(r, l) < 0)
                swap(l,r);
            }
    } else {
        if(compar(r,j) < 0){
            swap(j,r);
        } else if(compar(l, r) < 0)
            swap(l,r);
    }
    j = r ;
    for (;;) {
        while (compar(i+=size, v) < 0) ;
        while (compar(v, j-=size) < 0)
            if(j==i) break ;
        if(i>=j) break ;
        swap(i,j);
    }
    swap(i,r);
#pragma omp task
        omp_qsort (l, ((long long)i-(long long)l) / size, size, compar);
#pragma omp task
        omp_qsort (i+size, num - ((long long)i-(long long)l)/size - 1, size, compar);
}

void test_ego_loop3_long(size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, size_t *countresult, int stripes, int KBLOCK){
    size_t result = 0;

    int EGO_n = (n);
    int EGO_d = (d);
    double EGO_epsilon = (epsilon);
    double * EGO_array = (array);
    int EGO_blocks = (EGO_d + KBLOCK - 1) / KBLOCK;
    int EGO_stripes = (stripes);
    epsilonGridOrdering(EGO_n, EGO_d, EGO_epsilon, EGO_array, NUM_THREADS);
    int nn = ceilpowtwo(EGO_n);
    int **lower = (int **) malloc (EGO_stripes*sizeof(int*));
    int **upper = (int **) malloc (EGO_stripes*sizeof(int*));
    double *self = callocA64(sizeof (double) * EGO_n * EGO_blocks);
    prepareStripes(EGO_n, EGO_d, EGO_stripes, EGO_epsilon, EGO_array, lower, upper, (double *)0, NUM_THREADS);
    EGO_epsilon = EGO_epsilon * EGO_epsilon / 2;
    for (int i=0 ; i<EGO_n ; i++){
        double h=EGO_epsilon;
        for(int j=0 ; j<2*KBLOCK && j<EGO_d ; j++)
            h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];
        self[i/8*8*EGO_blocks+i%8]=h/2;
    }
    for(int k=2 ; k<EGO_blocks ; k++)
        for (int i=0 ; i<EGO_n ; i++){
            double h=0;
            for(int j=k * KBLOCK ; j<(k+1) * KBLOCK && j<EGO_d ; j++)
            h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];
        self[(i/8*EGO_blocks+k)*8+i%8]=h/2;
    }
    transpose_8xd(EGO_n, EGO_d, EGO_array,NUM_THREADS);
    // printf("end transpose8xd\n"); fflush(stdout);
    long long overall_load = 0;
    for (int i = 0; i < EGO_n / 8; i++)
        for (int j = 0; j < EGO_stripes; j++)
            overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];
    int loadstart[NUM_THREADS + 2];
    for (int i = 0; i <= NUM_THREADS; i++)
        loadstart[i] = 0;
    long long cum_load = 0;
    for (int i = 0; i < EGO_n / 8; i++) {
        for (int j = 0; j < EGO_stripes; j++)
            cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];
        loadstart[cum_load * NUM_THREADS / overall_load + 1] = i;
    }
    loadstart[NUM_THREADS] = EGO_n / 8;
    for (int i = 1; i <= NUM_THREADS; i++){
        if (loadstart[i] == 0)
            loadstart[i] = loadstart[i - 1];
    }

#pragma omp parallel for reduction(+:result)
    for (int par = 0; par < NUM_THREADS; par++) {
        veci64 resultvec = _mm512_setzero_si512();
        veci64 eights = _mm512_set1_epi64(8ll) ;

        int imin = loadstart[par];
        int imax = loadstart[par + 1];
        printf("STRIPES loop - %2d\n", omp_get_thread_num()); fflush(stdout);
        for (int s = 0; s < EGO_stripes; s++) {
            int i = 0;
            int j = 0;
            register veci64 const1 = _mm512_set1_epi64(1LL) ;
            register veci64 const2 = _mm512_add_epi64(const1, const1) ;
            register veci64 const3 = _mm512_add_epi64(const1, const2) ;
            register veci64 const4 = _mm512_add_epi64(const1, const3) ;
            register veci64 const5 = _mm512_add_epi64(const1, const4) ;
            register veci64 const6 = _mm512_add_epi64(const1, const5) ;
            register veci64 const7 = _mm512_add_epi64(const1, const6) ;
            register veci64 const0 = _mm512_setzero_si512();
            // Hilbert loop
            printf("before hilbert loop - %2d\n", omp_get_thread_num()); fflush(stdout);
            FGF_HILBERT_FOR(i, j, EGO_n / 8, EGO_n / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],
            FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&
            FURHIL_ub1 - 1 >= lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&
            FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {
                printf("in hilbert loop - %2d\n", omp_get_thread_num()); fflush(stdout);
                    register vec vi = _mm512_load_pd(self + i * 8 * EGO_blocks);
                    register vec vj = _mm512_load_pd(self + j * 8 * EGO_blocks);
                    register vec sum1 = vi + _mm512_permutexvar_pd(const0, vj);
                    register vec sum2 = vi + _mm512_permutexvar_pd(const1, vj);
                    register vec sum3 = vi + _mm512_permutexvar_pd(const2, vj);
                    register vec sum4 = vi + _mm512_permutexvar_pd(const3, vj);
                    register vec sum5 = vi + _mm512_permutexvar_pd(const4, vj);
                    register vec sum6 = vi + _mm512_permutexvar_pd(const5, vj);
                    register vec sum7 = vi + _mm512_permutexvar_pd(const6, vj);
                    register vec sum8 = vi + _mm512_permutexvar_pd(const7, vj);
                    vi = _mm512_load_pd(EGO_array + (EGO_d * i) * 8);
                    vj = _mm512_load_pd(EGO_array + (EGO_d * j) * 8);
                    sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);
                    sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);
                    sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);
                    sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);
                    sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);
                    sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);
                    sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);
                    sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);
                    int k; for (k=1 ; k<d ; k++){
                        if(k % KBLOCK == 0 && k > KBLOCK){
                            register veci64 allind = _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63);
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63);
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63);
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63);
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63);
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63);
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63);
                            allind += _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63);
                            if(_mm512_reduce_add_epi64(allind) >= 64) {k=d+1; break;}
                            vi = _mm512_load_pd(self + (i * EGO_blocks + k/KBLOCK) * 8);
                            vj = _mm512_load_pd(self + (j * EGO_blocks + k/KBLOCK) * 8);
                            sum1 += vi + _mm512_permutexvar_pd(const0, vj);
                            sum2 += vi + _mm512_permutexvar_pd(const1, vj);
                            sum3 += vi + _mm512_permutexvar_pd(const2, vj);
                            sum4 += vi + _mm512_permutexvar_pd(const3, vj);
                            sum5 += vi + _mm512_permutexvar_pd(const4, vj);
                            sum6 += vi + _mm512_permutexvar_pd(const5, vj);
                            sum7 += vi + _mm512_permutexvar_pd(const6, vj);
                            sum8 += vi + _mm512_permutexvar_pd(const7, vj);
                        }
                        vi = _mm512_load_pd(EGO_array + (EGO_d * i + k) * 8);
                        vj = _mm512_load_pd(EGO_array + (EGO_d * j + k) * 8);
                        sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);
                        sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);
                        sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);
                        sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);
                        sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);
                        sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);
                        sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);
                        sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);
                    }
                    if(k<=d){{
                        if (i==j){
                            sum1 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum1), 255, 0xbff0000000000000ull));
                            sum2 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum2), 254, 0xbff0000000000000ull));
                            sum3 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum3), 252, 0xbff0000000000000ull));
                            sum4 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum4), 248, 0xbff0000000000000ull));
                            sum5 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum5), 240, 0xbff0000000000000ull));
                            sum6 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum6), 224, 0xbff0000000000000ull));
                            sum7 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum7), 192, 0xbff0000000000000ull));
                            sum8 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum8), 128, 0xbff0000000000000ull));

                            resultvec += eights - _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63)         - _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63)
                                        - _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63)
                                        - _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63)
                                        - _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63)
                                        - _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63)
                                        - _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63)
                                        - _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63) ;
                        }
                }}
                FGF_HILBERT_END(i, j);
            }



            result += _mm512_reduce_add_epi64(resultvec);
//        double testres[8] __attribute__((aligned(64)));
//        _mm512_store_epi64(testres, resultvec);
//        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
        }

    }
    transpose_dx8(EGO_n, EGO_d, EGO_array,NUM_THREADS);

    // printf("result %zu\n", result);
    *countresult = result;
}

void test_ego_loop3_macro(size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, size_t *countresult, int stripes, int KBLOCK){
    size_t result = 0;

    EGO_PARALLEL_TRAN(n, d, epsilon, 5, array,NUM_THREADS)
        #pragma omp parallel for reduction(+:result)
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
        result += _mm512_reduce_add_epi64(resultvec);
//        double testres[8] __attribute__((aligned(64)));
//        _mm512_store_epi64(testres, resultvec);
//        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
    }
    EGO_END_TRAN

    *countresult = result;
}
