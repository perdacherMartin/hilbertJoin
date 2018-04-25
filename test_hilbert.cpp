void test_ego_loop3_long((size_t n, size_t d, size_t NUM_THREADS, double epsilon, double *array, size_t *countresult){
    size_t result = 0;

    int EGO_n = (n);
    int EGO_d = (d);
    double EGO_epsilon = (epsilon);
    double * EGO_array = (array);
    int EGO_blocks = (EGO_d + KBLOCK - 1) / KBLOCK;
    int EGO_stripes = (stripes);
    epsilonGridOrdering(EGO_n, EGO_d, EGO_epsilon, EGO_array);
    int nn = ceilpowtwo(EGO_n);
    int **lower = (int **) malloc (EGO_stripes*sizeof(int*));
    int **upper = (int **) malloc (EGO_stripes*sizeof(int*));
    double *self = callocA64(sizeof (double) * EGO_n * EGO_blocks);
    prepareStripes(EGO_n, EGO_d, EGO_stripes, EGO_epsilon, EGO_array, lower, upper, (double *)0);
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
    transpose_8xd(EGO_n, EGO_d, EGO_array);
    printf("end transpose8xd\n");
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
            FGF_HILBERT_FOR(i, j, EGO_n / 8, EGO_n / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],
            FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&
            FURHIL_ub1 - 1 >= lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&
            FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {
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
    transpose_dx8(EGO_n, EGO_d, EGO_array);

    // printf("result %zu\n", result);
    *countresult = result;
}
