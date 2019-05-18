
#ifndef HILLOOP_H
#define	HILLOOP_H

extern const long long HILLOOP_programList[][4];

#define FGF_HILBERT_FOR(i, j, imax, jmax, cond, cond2) {\
    register union {\
        double HILLOOP_hflo;\
        unsigned long long HILLOOP_hint;\
    };\
    int HILLOOP_imax = (imax);\
    int HILLOOP_jmax = (jmax);\
    unsigned long long HILLOOP_hilbert = 0ull;\
    int HILLOOP_c = 2;\
    int HILLOOP_maxext = imax - (i) > jmax - (j) ? imax - (i) : jmax-(j) ;\
    int HILLOOP_level;\
    for (HILLOOP_level = 0; HILLOOP_level < 32; HILLOOP_level++)\
        if (HILLOOP_maxext <= 1 << HILLOOP_level)\
            break;\
    int HILLOOP_maxvals = 1 <<HILLOOP_level ;\
    HILLOOP_level -= 3;\
    unsigned long long HILLOOP_maxhil = 1ULL <<(HILLOOP_level * 2);\
    while (HILLOOP_hilbert < HILLOOP_maxhil) {\
        int FURHIL_clevel ;\
        int FURHIL_lb0 = (i) ;\
        int FURHIL_ub0 = (i) + 1 ;\
        int FURHIL_lb1 = (j) ;\
        int FURHIL_ub1 = (j) + 1 ;\
        if (HILLOOP_c & 2){\
            FURHIL_ub0 += 3;\
            FURHIL_ub1 += 3;\
            for(FURHIL_clevel = 2 ; FURHIL_clevel < HILLOOP_level+3 ; FURHIL_clevel++){\
                FURHIL_ub0 += (1ULL<<FURHIL_clevel) ;\
                FURHIL_ub1 += (1ULL<<FURHIL_clevel) ;\
                if (FURHIL_lb0 < HILLOOP_imax && FURHIL_lb1 < HILLOOP_jmax && (cond2))\
                    break ;\
            }\
        } else {\
            FURHIL_lb0 -= 3;\
            FURHIL_lb1 -= 3;\
            for(FURHIL_clevel = 2 ; FURHIL_clevel < HILLOOP_level+3 ; FURHIL_clevel++){\
                FURHIL_lb0 -= (1ULL<<FURHIL_clevel) ;\
                FURHIL_lb1 -= (1ULL<<FURHIL_clevel) ;\
                if (FURHIL_lb0 < HILLOOP_imax && FURHIL_lb1 < HILLOOP_jmax && (cond2))\
                    break ;\
        }   }\
        if (FURHIL_clevel > 2) {\
            HILLOOP_c ^= (~FURHIL_clevel) & 1;\
            (j) += ((HILLOOP_c - 1) % 2) * ((1ULL << FURHIL_clevel) - 1);\
            (i) += ((HILLOOP_c - 2) % 2) * ((1ULL << FURHIL_clevel) - 1);\
            HILLOOP_hilbert += (1ULL << ((FURHIL_clevel-3) * 2)) ;\
            HILLOOP_c ^= 3*((~FURHIL_clevel) & 1);\
        } else {\
            unsigned long long HILLOOP_nanoH = HILLOOP_nanoprog[8][8][HILLOOP_c][0] ;\
            unsigned long long HILLOOP_nanoL = HILLOOP_nanoprog[8][8][HILLOOP_c][1] ;\
            for(;;){\
                if ((i)<HILLOOP_imax && (j)<HILLOOP_jmax && (cond))

#define FGF_HILBERT_END(i, j)\
                if(HILLOOP_nanoL == 1)\
                    break;\
                int HILLOOP_d = HILLOOP_nanoL & 1 | HILLOOP_nanoH & 2 ;\
                (i) += (HILLOOP_d - 2) % 2 ;\
                (j) += (HILLOOP_d - 1) % 2 ;\
                HILLOOP_nanoL /= 2 ;\
                HILLOOP_nanoH /= 2 ;\
            }\
            HILLOOP_hilbert++;\
        }\
        HILLOOP_hflo = HILLOOP_hilbert & (-HILLOOP_hilbert);\
        HILLOOP_level = (HILLOOP_hint - 4607182418800017408LL) >> 53;\
        int HILLOOP_test = HILLOOP_level & 1;\
        int HILLOOP_action = (HILLOOP_hilbert >> (2 * HILLOOP_level)) & 3;\
        HILLOOP_c ^= 3 * (HILLOOP_test != (HILLOOP_action == 3));\
        (j) += (HILLOOP_c - 1) % 2;\
        (i) += (HILLOOP_c - 2) % 2;\
        HILLOOP_c ^= (HILLOOP_test != (HILLOOP_action == 1));\
    }\
}

#define FGFN_HILBERT_FOR(i, j, imax, jmax, cond, cond2) {\
    register union {\
        double HILLOOP_hflo;\
        unsigned long long HILLOOP_hint;\
    };\
    int HILLOOP_imax = (imax);\
    int HILLOOP_jmax = (jmax);\
    unsigned long long HILLOOP_hilbert = 0ull;\
    int HILLOOP_c = 2; int HILLOOP_d;\
    int HILLOOP_maxext = imax - (i) > jmax - (j) ? imax - (i) : jmax-(j) ;\
    int HILLOOP_level;\
    for (HILLOOP_level = 0; HILLOOP_level < 32; HILLOOP_level++)\
        if (HILLOOP_maxext <= 1 << HILLOOP_level)\
            break;\
    int HILLOOP_maxvals = 1<<HILLOOP_level ;\
    HILLOOP_level -= 3;\
    unsigned long long HILLOOP_maxhil = 1<<(HILLOOP_level * 2);\
    while (HILLOOP_hilbert < HILLOOP_maxhil) {\
        int FURHIL_clevel ;\
        int FURHIL_lb0 = (i) ;\
        int FURHIL_ub0 = (i) + 1 ;\
        int FURHIL_lb1 = (j) ;\
        int FURHIL_ub1 = (j) + 1 ;\
        if (HILLOOP_c & 2){\
            FURHIL_ub0 += 3;\
            FURHIL_ub1 += 3;\
            for(FURHIL_clevel = 2 ; FURHIL_clevel < HILLOOP_level+3 ; FURHIL_clevel++){\
                FURHIL_ub0 += (1<<FURHIL_clevel) ;\
                FURHIL_ub1 += (1<<FURHIL_clevel) ;\
                if (FURHIL_lb0 < HILLOOP_imax && FURHIL_lb1 < HILLOOP_jmax && (cond2))\
                    break ;\
            }\
        } else {\
            FURHIL_lb0 -= 3;\
            FURHIL_lb1 -= 3;\
            for(FURHIL_clevel = 2 ; FURHIL_clevel < HILLOOP_level+3 ; FURHIL_clevel++){\
                FURHIL_lb0 -= (1<<FURHIL_clevel) ;\
                FURHIL_lb1 -= (1<<FURHIL_clevel) ;\
                if (FURHIL_lb0 < HILLOOP_imax && FURHIL_lb1 < HILLOOP_jmax && (cond2))\
                    break ;\
        }   }\
        if (FURHIL_clevel > 2) {\
            HILLOOP_c ^= (~FURHIL_clevel) & 1;\
            (j) += ((HILLOOP_c - 1) % 2) * ((1 << FURHIL_clevel) - 1);\
            (i) += ((HILLOOP_c - 2) % 2) * ((1 << FURHIL_clevel) - 1);\
            HILLOOP_hilbert += (1 << ((FURHIL_clevel-3) * 2)) ;\
            HILLOOP_c ^= 3*((~FURHIL_clevel) & 1);\
        } else {\
            unsigned long long HILLOOP_nano = HILLOOP_fgfnano[HILLOOP_c][0] ;\
            for(;;){\
                if ((i)<HILLOOP_imax && (j)<HILLOOP_jmax && (cond))

#define FGFN_HILBERT_END(i, j)\
                if(HILLOOP_nano < 3){\
                    if(HILLOOP_nano == 1)\
                        break;\
                    HILLOOP_nano = HILLOOP_fgfnano[HILLOOP_c][1];\
                    HILLOOP_d = HILLOOP_c;\
                }else{\
                    HILLOOP_d = HILLOOP_nano & 3 ;\
                    HILLOOP_nano /= 4 ;\
                }\
                (i) += (HILLOOP_d - 2) % 2 ;\
                (j) += (HILLOOP_d - 1) % 2 ;\
            }\
            HILLOOP_hilbert++;\
        }\
        HILLOOP_hflo = HILLOOP_hilbert & (-HILLOOP_hilbert);\
        HILLOOP_level = (HILLOOP_hint - 4607182418800017408LL) >> 53;\
        int HILLOOP_test = HILLOOP_level & 1;\
        int HILLOOP_action = (HILLOOP_hilbert >> (2 * HILLOOP_level)) & 3;\
        HILLOOP_c ^= 3 * (HILLOOP_test != (HILLOOP_action == 3));\
        (j) += (HILLOOP_c - 1) % 2;\
        (i) += (HILLOOP_c - 2) % 2;\
        HILLOOP_c ^= (HILLOOP_test != (HILLOOP_action == 1));\
    }\
}

#define FGFC_HILBERT_FOR(i, j, imax, jmax, cond, cond2, cond3) {\
    register union {\
        double HILLOOP_hflo;\
        unsigned long long HILLOOP_hint;\
    };\
    int HILLOOP_imax = (imax);\
    int HILLOOP_jmax = (jmax);\
    unsigned long long HILLOOP_hilbert = 0ull;\
    int HILLOOP_c = 2; int HILLOOP_d;\
    int HILLOOP_maxext = imax - (i) > jmax - (j) ? imax - (i) : jmax-(j) ;\
    int HILLOOP_level;\
    for (HILLOOP_level = 0; HILLOOP_level < 32; HILLOOP_level++)\
        if (HILLOOP_maxext <= 1 << HILLOOP_level)\
            break;\
    int HILLOOP_maxvals = 1<<HILLOOP_level ;\
    HILLOOP_level -= 3;\
    unsigned long long HILLOOP_maxhil = 1<<(HILLOOP_level * 2);\
    while (HILLOOP_hilbert < HILLOOP_maxhil) {\
        int FURHIL_clevel ;\
        int FURHIL_lb0 = (i) ;\
        int FURHIL_ub0 = (i) + 1 ;\
        int FURHIL_lb1 = (j) ;\
        int FURHIL_ub1 = (j) + 1 ;\
        if (HILLOOP_c & 2){\
            FURHIL_ub0 += 3;\
            FURHIL_ub1 += 3;\
            for(FURHIL_clevel = 2 ; FURHIL_clevel < HILLOOP_level+3 ; FURHIL_clevel++){\
                FURHIL_ub0 += (1<<FURHIL_clevel) ;\
                FURHIL_ub1 += (1<<FURHIL_clevel) ;\
                if (FURHIL_lb0 < HILLOOP_imax && FURHIL_lb1 < HILLOOP_jmax && (cond2))\
                    break ;\
            }\
        } else {\
            FURHIL_lb0 -= 3;\
            FURHIL_lb1 -= 3;\
            for(FURHIL_clevel = 2 ; FURHIL_clevel < HILLOOP_level+3 ; FURHIL_clevel++){\
                FURHIL_lb0 -= (1<<FURHIL_clevel) ;\
                FURHIL_lb1 -= (1<<FURHIL_clevel) ;\
                if (FURHIL_lb0 < HILLOOP_imax && FURHIL_lb1 < HILLOOP_jmax && (cond2))\
                    break ;\
        }   }\
        if (FURHIL_clevel > 2) {\
            HILLOOP_c ^= (~FURHIL_clevel) & 1;\
            (j) += ((HILLOOP_c - 1) % 2) * ((1 << FURHIL_clevel) - 1);\
            (i) += ((HILLOOP_c - 2) % 2) * ((1 << FURHIL_clevel) - 1);\
            HILLOOP_hilbert += (1 << ((FURHIL_clevel-3) * 2)) ;\
            HILLOOP_c ^= 3*((~FURHIL_clevel) & 1);\
        } else {\
            unsigned long long HILLOOP_nano = HILLOOP_fgfnano[HILLOOP_c][0] ;\
            FURHIL_clevel = 2;\
            FURHIL_lb0 = (i) ;\
            FURHIL_ub0 = (i) + 1 ;\
            FURHIL_lb1 = (j) ;\
            FURHIL_ub1 = (j) + 1 ;\
            if (HILLOOP_c & 2){\
                FURHIL_ub0 += 7;\
                FURHIL_ub1 += 7;\
            } else {\
                FURHIL_lb0 -= 7;\
                FURHIL_lb1 -= 7;\
            }\
            if(FURHIL_ub0 >= HILLOOP_imax || FURHIL_ub1 >= HILLOOP_jmax || !(cond3)){\
                for(;;){\
                    if ((i)<HILLOOP_imax && (j)<HILLOOP_jmax && (cond))

#define FGFC_HILBERT_ELSE(i, j)\
                    if(HILLOOP_nano < 3){\
                        if(HILLOOP_nano == 1)\
                            break;\
                        HILLOOP_nano = HILLOOP_fgfnano[HILLOOP_c][1];\
                        HILLOOP_d = HILLOOP_c;\
                    }else{\
                        HILLOOP_d = HILLOOP_nano & 3 ;\
                        HILLOOP_nano /= 4 ;\
                    }\
                    (i) += (HILLOOP_d - 2) % 2 ;\
                    (j) += (HILLOOP_d - 1) % 2 ;\
                }\
                HILLOOP_hilbert++;\
            } else {\
                for (;;){

#define FGFC_HILBERT_END(i, j)\
                    if(HILLOOP_nano < 3){\
                        if(HILLOOP_nano == 1)\
                            break;\
                        HILLOOP_nano = HILLOOP_fgfnano[HILLOOP_c][1];\
                        HILLOOP_d = HILLOOP_c;\
                    }else{\
                        HILLOOP_d = HILLOOP_nano & 3 ;\
                        HILLOOP_nano /= 4 ;\
                    }\
                    (i) += (HILLOOP_d - 2) % 2 ;\
                    (j) += (HILLOOP_d - 1) % 2 ;\
                }\
                HILLOOP_hilbert++;\
            }\
        }\
        HILLOOP_hflo = HILLOOP_hilbert & (-HILLOOP_hilbert);\
        HILLOOP_level = (HILLOOP_hint - 4607182418800017408LL) >> 53;\
        int HILLOOP_test = HILLOOP_level & 1;\
        int HILLOOP_action = (HILLOOP_hilbert >> (2 * HILLOOP_level)) & 3;\
        HILLOOP_c ^= 3 * (HILLOOP_test != (HILLOOP_action == 3));\
        (j) += (HILLOOP_c - 1) % 2;\
        (i) += (HILLOOP_c - 2) % 2;\
        HILLOOP_c ^= (HILLOOP_test != (HILLOOP_action == 1));\
    }\
}

//                if(HILLOOP_nano < 3){\
//                    if(HILLOOP_nano == 1)\
//                        break;\
//                    HILLOOP_nano = HILLOOP_fgfnano[HILLOOP_c][1];\
//                    (i) += (HILLOOP_c - 2) % 2 ;\
//                    (j) += (HILLOOP_c - 1) % 2 ;\
//                    continue;\
//                }\
//                HILLOOP_d = HILLOOP_nano & 3 ;\
//                HILLOOP_nano /= 4 ;\
//                (i) += (HILLOOP_d - 2) % 2 ;\
//                (j) += (HILLOOP_d - 1) % 2 ;\
//            }\


#endif	/* HILLOOP_H */
