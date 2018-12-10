
#ifndef HILLOOP_H
#define	HILLOOP_H

extern const long long HILLOOP_programList[][4];

#define HILLOOP_START(i, j, imin, imax, jmin, jmax)\
            if ((imin) < (imax) && (jmin) < (jmax)) {\
                int HILLOOP_imin = (imin) ;\
                int HILLOOP_imax = (imax) ;\
                int HILLOOP_jmin = (jmin) ;\
                int HILLOOP_jmax = (jmax) ;\
                (i) = HILLOOP_imin;\
                (j) = HILLOOP_jmin;\
                int HILLOOP_c;\
                unsigned long long HILLOOP_hilbert = 0;\
                register union {\
                        double HILLOOP_hflo;\
                        unsigned long long HILLOOP_hint;\
                };\
                unsigned long long HILLOOP_twopot, HILLOOP_tt;\
                int HILLOOP_in, HILLOOP_ix, HILLOOP_jn, HILLOOP_jx, HILLOOP_creset;\
                long long HILLOOP_program;\
                if (HILLOOP_imax - HILLOOP_imin < HILLOOP_jmax - HILLOOP_jmin) {\
                    if (HILLOOP_imax - HILLOOP_imin == 1) {\
                        HILLOOP_tt = 2;\
                        HILLOOP_twopot = 1;\
                        HILLOOP_creset = HILLOOP_c = 2;\
                    } else {\
                        HILLOOP_hflo = HILLOOP_imax - HILLOOP_imin;\
                        HILLOOP_hint = (HILLOOP_hint >> 52) - 1023;\
                        HILLOOP_tt = 1 << HILLOOP_hint;\
                        HILLOOP_twopot = HILLOOP_tt / 2;\
                        HILLOOP_creset = HILLOOP_c = 3 - (HILLOOP_hint & 1);\
                    }\
                    HILLOOP_jn = HILLOOP_tt + (HILLOOP_jmax - HILLOOP_jmin) % HILLOOP_tt < HILLOOP_jmax - HILLOOP_jmin ? HILLOOP_tt + (HILLOOP_jmax - HILLOOP_jmin) % HILLOOP_tt : HILLOOP_jmax - HILLOOP_jmin;\
                    HILLOOP_ix = HILLOOP_imax;\
                    HILLOOP_in = HILLOOP_ix - HILLOOP_imin;\
                    HILLOOP_jx = HILLOOP_jn + HILLOOP_jmin;\
                } else {\
                    if (HILLOOP_jmax - HILLOOP_jmin == 1) {\
                        HILLOOP_tt = 2;\
                        HILLOOP_twopot = 1;\
                        HILLOOP_creset = HILLOOP_c = 3;\
                    } else {\
                        HILLOOP_hflo = HILLOOP_jmax - HILLOOP_jmin;\
                        HILLOOP_hint = (HILLOOP_hint >> 52) - 1023;\
                        HILLOOP_tt = 1 << HILLOOP_hint;\
                        HILLOOP_twopot = HILLOOP_tt / 2;\
                        HILLOOP_creset = HILLOOP_c = 2 + (HILLOOP_hint & 1);\
                    }\
                    HILLOOP_in = HILLOOP_tt + (HILLOOP_imax-HILLOOP_imin) % HILLOOP_tt < HILLOOP_imax-HILLOOP_imin ? HILLOOP_tt + (HILLOOP_imax-HILLOOP_imin) % HILLOOP_tt : HILLOOP_imax-HILLOOP_imin;\
                    HILLOOP_jx = HILLOOP_jmax;\
                    HILLOOP_jn = HILLOOP_jx - HILLOOP_jmin;\
                    HILLOOP_ix = HILLOOP_in + HILLOOP_imin;\
                }\
                HILLOOP_program = HILLOOP_programList[((HILLOOP_in + HILLOOP_twopot - 1) / HILLOOP_twopot)*5 + ((HILLOOP_jn + HILLOOP_twopot - 1) / HILLOOP_twopot)][HILLOOP_c];\
                while (1) {

#define HILLOOP_END(i,j)\
                    if (HILLOOP_program) {\
                        (i) += (HILLOOP_program & 3) - 2;\
                        HILLOOP_program >>= 2;\
                        (j) += (HILLOOP_program & 3) - 2;\
                        HILLOOP_program >>= 2;\
                    } else {\
                        HILLOOP_hilbert++;\
                        HILLOOP_hflo = HILLOOP_hilbert & -HILLOOP_hilbert;\
                        int HILLOOP_test = (HILLOOP_hint - 4607182418800017408LL) >> 53; \
                        int HILLOOP_action = (HILLOOP_hilbert >> (2 * HILLOOP_test)) & 3;\
                        HILLOOP_test &= 1;\
                        HILLOOP_c ^= 3 * (HILLOOP_test != (HILLOOP_action == 3));\
                        (j) += (HILLOOP_c - 1) % 2;\
                        (i) -= (2 - HILLOOP_c) % 2;\
                        HILLOOP_c ^= (HILLOOP_test != (HILLOOP_action == 1));\
                        if ((i) == HILLOOP_ix) {\
                            if ((i) >= HILLOOP_imax)\
                                break;\
                            HILLOOP_in = HILLOOP_tt;\
                            HILLOOP_ix += HILLOOP_tt;\
                            HILLOOP_hilbert = 0;\
                            HILLOOP_c = HILLOOP_creset;\
                        }\
                        if ((j) == HILLOOP_jx) {\
                            if ((j) >= HILLOOP_jmax)\
                                break;\
                            HILLOOP_jn = HILLOOP_tt;\
                            HILLOOP_jx += HILLOOP_tt;\
                            HILLOOP_hilbert = 0;\
                            HILLOOP_c = HILLOOP_creset;\
                        }\
                        HILLOOP_program = HILLOOP_programList[5 * ((((((i)-HILLOOP_imin) * HILLOOP_twopot / HILLOOP_in) + 1) * HILLOOP_in + HILLOOP_twopot - 1) / HILLOOP_twopot - ((((i)-HILLOOP_imin) * HILLOOP_twopot / HILLOOP_in) * HILLOOP_in + HILLOOP_twopot - 1) / HILLOOP_twopot)\
                                +((((((j)-HILLOOP_jmin) * HILLOOP_twopot / HILLOOP_jn) + 1) * HILLOOP_jn + HILLOOP_twopot - 1) / HILLOOP_twopot - ((((j)-HILLOOP_jmin) * HILLOOP_twopot / HILLOOP_jn) * HILLOOP_jn + HILLOOP_twopot - 1) / HILLOOP_twopot)][HILLOOP_c];\
                    }\
                }\
            }

#define CANOLOOP_START(i, j, imin, imax, jmin, jmax) {\
            for ((i)=(imin) ; (i)<(imax) ; (i)++)\
                for ((j)=(jmin) ; (j)<(jmax) ; (j)++) {

#define CANOLOOP_END(i,j) } }

#define CANO2LOOP_START(i, j, imin, imax, jmin, jmax) {\
            for ((j)=(jmin) ; (j)<(jmax) ; (j)++)\
                for ((i)=(imin) ; (i)<(imax) ; (i)++){

#define CANO2LOOP_END(i,j) } }

#define CCILOOP_START(i,j,imin,imax,jmin,jmax,istep) {\
            for (int HILLOOP_I=(imin) ; HILLOOP_I<(imax) ; HILLOOP_I+=(istep))\
                for ((j)=(jmin) ; (j)<(jmax) ; (j)++)\
                    for ((i)=HILLOOP_I ; (i)<HILLOOP_I+(istep) && (i)<(imax) ; (i)++) {

#define CCILOOP_END(i,j) } }

#define CCJLOOP_START(i,j,imin,imax,jmin,jmax,jstep) {\
            for (int HILLOOP_J=(jmin) ; HILLOOP_J<(jmax) ; HILLOOP_J+=(jstep))\
                for ((i)=(imin) ; (i)<(imax) ; (i)++)\
                    for ((j)=HILLOOP_J ; (j)<HILLOOP_J+(jstep) && (j)<(jmax) ; (j)++) {

#define CCJLOOP_END(i,j) } }

#define CCIJLOOP_START(i,j,imin,imax,jmin,jmax,istep,jstep) {\
            for (int HILLOOP_I=(imin) ; HILLOOP_I<(imax) ; HILLOOP_I+=(istep))\
                for (int HILLOOP_J=(jmin) ; HILLOOP_J<(jmax) ; HILLOOP_J+=(jstep))\
                    for ((i)=HILLOOP_I ; (i)<HILLOOP_I+(istep) && (i)<(imax) ; (i)++)\
                        for ((j)=HILLOOP_J ; (j)<HILLOOP_J+(jstep) && (j)<(jmax) ; (j)++) {

#define CCIJLOOP_END(i,j) } }

extern const unsigned long long HILLOOP_mortonListI[];

extern const unsigned long long HILLOOP_mortonListJ[];

#define ZORDLOOP_START(i, j, imin, imax, jmin, jmax) {\
            int HILLOOP_imin = (imin) ;\
            int HILLOOP_jmin = (jmin) ;\
            int HILLOOP_idiff = (imax) - HILLOOP_imin;\
            int HILLOOP_jdiff = (jmax) - HILLOOP_jmin;\
            if (HILLOOP_idiff > 0 && HILLOOP_jdiff > 0) {\
                union {\
                    double HILLOOP_ddd;\
                    long long HILLOOP_lll;\
                };\
                HILLOOP_ddd = HILLOOP_idiff < HILLOOP_jdiff ? HILLOOP_idiff : HILLOOP_jdiff;\
                unsigned long long HILLOOP_t = (HILLOOP_lll >> 52) - 1023;\
                if (HILLOOP_t) HILLOOP_t--;\
                unsigned long long HILLOOP_tt = 1 << HILLOOP_t;\
                unsigned long long HILLOOP_tti = (HILLOOP_idiff / HILLOOP_tt / 2) * HILLOOP_tt ;\
                unsigned long long HILLOOP_ttj = (HILLOOP_jdiff / HILLOOP_tt / 2) * HILLOOP_tt ;\
                HILLOOP_tti = HILLOOP_tti<1 ? 1 : HILLOOP_tti ;\
                HILLOOP_ttj = HILLOOP_ttj<1 ? 1 : HILLOOP_ttj ;\
                unsigned long long HILLOOP_I = 0;\
                unsigned long long HILLOOP_J = 0;\
                for (int HILLOOP_outer = 0; HILLOOP_outer < (HILLOOP_tti > HILLOOP_ttj ? HILLOOP_tti : HILLOOP_ttj) ; HILLOOP_outer += HILLOOP_tt) {\
                    if (HILLOOP_idiff < HILLOOP_jdiff) {\
                        HILLOOP_I = 0;\
                        HILLOOP_J = HILLOOP_outer ;\
                    }else{\
                        HILLOOP_J = 0;\
                        HILLOOP_I = HILLOOP_outer ;\
                    }\
                    unsigned long long HILLOOP_Z = 0;\
                    while (HILLOOP_Z < HILLOOP_tt * HILLOOP_tt) {\
                        (i) = HILLOOP_I * HILLOOP_idiff / HILLOOP_tti;\
                        (j) = HILLOOP_J * HILLOOP_jdiff / HILLOOP_ttj;\
                        int HILLOOP_Idiff = (HILLOOP_I + 1) * HILLOOP_idiff / HILLOOP_tti - (i);\
                        int HILLOOP_Jdiff = (HILLOOP_J + 1) * HILLOOP_jdiff / HILLOOP_ttj - (j);\
                        (i) += HILLOOP_imin;\
                        (j) += HILLOOP_jmin;\
                        unsigned long long HILLOOP_programI = HILLOOP_mortonListI[5 * HILLOOP_Idiff + HILLOOP_Jdiff];\
                        unsigned long long HILLOOP_programJ = HILLOOP_mortonListJ[5 * HILLOOP_Idiff + HILLOOP_Jdiff];\
                        do {

#define ZORDLOOP_END(i,j)\
                            (i) = (i) + (HILLOOP_programI & 7) - 5;\
                            (j) = (j) + (HILLOOP_programJ & 7) - 5;\
                            HILLOOP_programI >>= 3;\
                            HILLOOP_programJ >>= 3;\
                        } while (HILLOOP_programI);\
                        HILLOOP_Z++;\
                        HILLOOP_ddd = HILLOOP_Z & -HILLOOP_Z;\
                        int HILLOOP_action = (HILLOOP_lll >> 52) - 1023;\
                        if (HILLOOP_action & 1) {\
                            HILLOOP_J &= (-1 << ((HILLOOP_action + 1) / 2));\
                            HILLOOP_I++;\
                        } else {\
                            HILLOOP_I &= (-1 << (HILLOOP_action / 2));\
                            HILLOOP_J++;\
                        }\
                    }\
                }\
            }\
        }

extern const unsigned long long HILLOOP_asuzListI[][4];

extern const unsigned long long HILLOOP_asuzListJ[][4];

#define ASUZLOOP_START(i, j, imin, imax, jmin, jmax) {\
            int HILLOOP_imin = (imin) ;\
            int HILLOOP_jmin = (jmin) ;\
            int HILLOOP_idiff = (imax) - HILLOOP_imin;\
            int HILLOOP_jdiff = (jmax) - HILLOOP_jmin;\
            if (HILLOOP_idiff > 0 && HILLOOP_jdiff > 0) {\
                union {\
                    double HILLOOP_ddd;\
                    long long HILLOOP_lll;\
                };\
                HILLOOP_ddd = HILLOOP_idiff < HILLOOP_jdiff ? HILLOOP_idiff : HILLOOP_jdiff;\
                unsigned long long HILLOOP_t = (HILLOOP_lll >> 52) - 1023;\
                if (HILLOOP_t) HILLOOP_t--;\
                int HILLOOP_creset = 2 + ((HILLOOP_t & 1) == (HILLOOP_idiff < HILLOOP_jdiff));\
                unsigned long long HILLOOP_tt = 1 << HILLOOP_t;\
                unsigned long long HILLOOP_tti = (HILLOOP_idiff / HILLOOP_tt / 2) * HILLOOP_tt ;\
                unsigned long long HILLOOP_ttj = (HILLOOP_jdiff / HILLOOP_tt / 2) * HILLOOP_tt ;\
                HILLOOP_tti = HILLOOP_tti<1 ? 1 : HILLOOP_tti ;\
                HILLOOP_ttj = HILLOOP_ttj<1 ? 1 : HILLOOP_ttj ;\
                unsigned long long HILLOOP_I = 0;\
                unsigned long long HILLOOP_J = 0;\
                (i) = HILLOOP_imin ;\
                for (int HILLOOP_outer = 0; HILLOOP_outer < (HILLOOP_tti > HILLOOP_ttj ? HILLOOP_tti : HILLOOP_ttj) ; HILLOOP_outer += HILLOOP_tt) {\
                    if (HILLOOP_idiff < HILLOOP_jdiff) {\
                        HILLOOP_I = 0;\
                        HILLOOP_J = HILLOOP_outer ;\
                    }else{\
                        HILLOOP_J = 0;\
                        HILLOOP_I = HILLOOP_outer ;\
                    }\
                    int HILLOOP_c = HILLOOP_creset ;\
                    unsigned long long HILLOOP_Z = 0;\
                    while (HILLOOP_Z < HILLOOP_tt * HILLOOP_tt) {\
                        int HILLOOP_Idiff = (HILLOOP_I + 1) * HILLOOP_idiff / HILLOOP_tti - HILLOOP_I * HILLOOP_idiff / HILLOOP_tti;\
                        int HILLOOP_Jdiff = (HILLOOP_J + 1) * HILLOOP_jdiff / HILLOOP_ttj - HILLOOP_J * HILLOOP_jdiff / HILLOOP_ttj;\
                        (j) = HILLOOP_J * HILLOOP_jdiff / HILLOOP_ttj + HILLOOP_jmin;\
                        unsigned long long HILLOOP_programI = HILLOOP_asuzListI[5 * HILLOOP_Idiff + HILLOOP_Jdiff][HILLOOP_c];\
                        unsigned long long HILLOOP_programJ = HILLOOP_asuzListJ[5 * HILLOOP_Idiff + HILLOOP_Jdiff][HILLOOP_c];\
                        do {

#define ASUZLOOP_END(i,j)\
                            (i) = (i) + (HILLOOP_programI & 7) - 5;\
                            (j) = (j) + (HILLOOP_programJ & 7) - 5;\
                            HILLOOP_programI >>= 3;\
                            HILLOOP_programJ >>= 3;\
                        } while (HILLOOP_programI);\
                        HILLOOP_Z++;\
                        HILLOOP_ddd = HILLOOP_Z & -HILLOOP_Z;\
                        int HILLOOP_test = (HILLOOP_lll - 4607182418800017408LL) >> 53; \
                        int HILLOOP_action = (HILLOOP_Z >> (2 * HILLOOP_test)) & 3;\
                        HILLOOP_c ^= 3 * ((HILLOOP_test&1) != (HILLOOP_action == 3));\
                        if(HILLOOP_action == 2) {if (HILLOOP_c & 1) HILLOOP_J &= -2<<HILLOOP_test; else HILLOOP_J ++ ;}\
                        else {if (HILLOOP_c & 1) HILLOOP_J &= -1<<HILLOOP_test ; else HILLOOP_J ++;}\
                        HILLOOP_I -= (2 - HILLOOP_c) % 2;\
                        (i) -= (2 - HILLOOP_c) % 2;\
                        HILLOOP_c ^= ((HILLOOP_test&1) != (HILLOOP_action == 1));\
                    }\
                }\
            }\
        }

// Action Table
// State Action Follow i     j      d0 XOR a1   d1 XOR a0    d1 AND a0
// 00    00     01     +1    +1     0           0
// 00    01     00     0     +1     0           1
// 00    10     00     0     0      1           0
// 00    11     11     +1    0      1           1
// 01    00     00     +1    +1     1           0
// 01    01     01     +1    0      1           1
// 01    10     01     0     0      0           0
// 01    11     10     0     +1     0           1
// 10    00     11     0     0      0           1
// 10    01     10     +1    0      0           0
// 10    10     10     +1    +1     1           1
// 10    11     01     0     +1     1           0
// 11    00     10     0     0      1           1
// 11    01     11     0     +1     1           0
// 11    10     11     +1    +1     0           1
// 11    11     00     +1    0      0           0

#define BIGHLOOP_START(i, j, imin, imax, jmin, jmax) {\
    register union {\
        double HILLOOP_hflo;\
        unsigned long long HILLOOP_hint;\
    };\
    HILLOOP_hflo = ((imax)-(imin)>(jmax)-(jmin) ? (imax)-(imin) : (jmax)-(jmin)) - 1 ;\
    HILLOOP_hint = ((HILLOOP_hint >> 52) - 1022) * 2;\
    unsigned long long HILLOOP_hmax = 1 << HILLOOP_hint ;\
    for (unsigned long long HILLOOP_hilbert = 0; HILLOOP_hilbert < HILLOOP_hmax ; HILLOOP_hilbert ++) {\
        (i) = 0 ;\
        (j) = 0 ;\
        int HILLOOP_d = 3 ;\
        for (int HILLOOP_kk=HILLOOP_hint-2 ; HILLOOP_kk >= 0 ; HILLOOP_kk-=2) {\
            int HILLOOP_action = (HILLOOP_hilbert >> HILLOOP_kk) & 3 ;\
            (i) = 2 * (i) + (((HILLOOP_d | ~HILLOOP_action) & 1) ^ ((HILLOOP_d ^ HILLOOP_action)>>1)) ;\
            (j) = 2 * (j) + ((1 & ~(HILLOOP_d & HILLOOP_action)) ^ ((HILLOOP_d ^ HILLOOP_action)>>1)) ;\
            HILLOOP_d ^= (HILLOOP_action == 0);\
            HILLOOP_d ^= 3*(HILLOOP_action == 3) ;\
        }\
        (i) += (imin) ;\
        (j) += (jmin) ;\
        if ((i)<(imax) && (j)<(jmax)) {\

#define BIGHLOOP_END(i,j)\
        }\
    }\
}\

//int HILLOOP_peano_i[][9] = {{0,1,2,2,1,0,0,1,2},{0,1,2,2,1,0,0,1,2},{2,1,0,0,1,2,2,1,0},{2,1,0,0,1,2,2,1,0}} ;
//int HILLOOP_peano_j[][9] = {{0,0,0,1,1,1,2,2,2},{2,2,2,1,1,1,0,0,0},{0,0,0,1,1,1,2,2,2},{2,2,2,1,1,1,0,0,0}} ;
//int HILLOOP_peano_d[][9] = {{0,1,0,2,3,2,0,1,0},{1,0,1,3,2,3,1,0,1},{2,3,2,0,1,0,2,3,2},{3,2,3,1,0,1,3,2,3}} ;
//
//#define PEANOLOOP_START(i, j, imin, imax, jmin, jmax) {\
//    unsigned long long HILLOOP_pow = 1;\
//    while (HILLOOP_pow<imax-imin || HILLOOP_pow<jmax-jmin) HILLOOP_pow *= 3;\
//    HILLOOP_pow *= HILLOOP_pow;\
//    for (unsigned long long HILLOOP_hilbert=0 ; HILLOOP_hilbert<HILLOOP_pow ; HILLOOP_hilbert++) {\
//        (i) = 0;\
//        (j) = 0;\
//        int HILLOOP_d = 0 ;\
//        for (int HILLOOP_kk = HILLOOP_pow/9 ; HILLOOP_kk ; HILLOOP_kk/=9){\
//            int HILLOOP_a = (HILLOOP_hilbert/HILLOOP_kk) % 9 ;\
//            (i) = 3 * (i) + HILLOOP_peano_i[HILLOOP_d][HILLOOP_a] ;\
//            (j) = 3 * (j) + HILLOOP_peano_j[HILLOOP_d][HILLOOP_a] ;\
//            HILLOOP_d = HILLOOP_peano_d[HILLOOP_d][HILLOOP_a] ;\
//        }\
//        (i) += (imin) ;\
//        (j) += (jmin) ;\
//        if ((i)<(imax) && (j)<(jmax)) {\
//
//# define PEANOLOOP_END(i,j)\
//        }\
//    }\
//}\

#define PEANOLOOP_START(i, j, imin, imax, jmin, jmax) {\
    unsigned long long HILLOOP_pow = 1;\
    while (HILLOOP_pow<imax-imin || HILLOOP_pow<jmax-jmin) HILLOOP_pow *= 3;\
    HILLOOP_pow *= HILLOOP_pow;\
    (i) = (imin);\
    (j) = (jmin);\
    int HILLOOP_lr=1; int HILLOOP_ud=1;\
    for (unsigned long long HILLOOP_hilbert=1 ; HILLOOP_hilbert<=HILLOOP_pow ; HILLOOP_hilbert++) {\
        if((i)<(imax) && (j)<(jmax)) {

#define PEANOLOOP_END(i,j)\
        }\
        int HILLOOP_a = HILLOOP_hilbert ;\
        while(HILLOOP_a % 9 == 0) {\
            HILLOOP_a /= 9 ;\
        }\
        if(HILLOOP_a % 3){\
            (j) += HILLOOP_lr;\
            HILLOOP_ud = (-HILLOOP_ud);\
        } else {\
            (i) += HILLOOP_ud;\
            HILLOOP_lr = (-HILLOOP_lr);\
        }\
    }\
}

#define BIGZLOOP_START(i,j,imin,imax,jmin,jmax) {\
    register union {\
        double HILLOOP_hflo;\
        unsigned long long HILLOOP_hint;\
    };\
    HILLOOP_hflo = ((imax)-(imin)>(jmax)-(jmin) ? (imax)-(imin) : (jmax)-(jmin)) - 1 ;\
    HILLOOP_hint = ((HILLOOP_hint >> 52) - 1022) * 2;\
    unsigned long long HILLOOP_hmax = 1 << HILLOOP_hint ;\
    for (unsigned long long HILLOOP_hilbert = 0; HILLOOP_hilbert < HILLOOP_hmax ; HILLOOP_hilbert ++) {\
        (i) = 0 ;\
        (j) = 0 ;\
        for (int HILLOOP_kk=HILLOOP_hint-1 ; HILLOOP_kk >= 0 ; HILLOOP_kk--) {\
            (i) = 2 * (i) + ((HILLOOP_hilbert >> HILLOOP_kk) & 1) ;\
            HILLOOP_kk-- ;\
            (j) = 2 * (j) + ((HILLOOP_hilbert >> HILLOOP_kk) & 1) ;\
        }\
        (i) += (imin) ;\
        (j) += (jmin) ;\
        if ((i)<(imax) && (j)<(jmax)) {\

#define BIGZLOOP_END(i,j)\
        }\
    }\
}\

extern const unsigned long long HILLOOP_nanoprog[9][9][4][2];

extern const unsigned long long HILLOOP_fgfnano[4][2];

#define FUR_HILBERT_START(i,j,imin,imax,jmin,jmax) {\
    int HILLOOP_imin = (imin);\
    int HILLOOP_imax = (imax);\
    int HILLOOP_jmin = (jmin);\
    int HILLOOP_jmax = (jmax);\
    (i) = HILLOOP_imin;\
    (j) = HILLOOP_jmin;\
    int HILLOOP_idiff = HILLOOP_imax - HILLOOP_imin;\
    int HILLOOP_jdiff = HILLOOP_jmax - HILLOOP_jmin;\
    if (HILLOOP_idiff > 0 && HILLOOP_jdiff > 0) {\
        int HILLOOP_iceil = HILLOOP_idiff - 1;\
        HILLOOP_iceil |= HILLOOP_iceil >> 1;\
        HILLOOP_iceil |= HILLOOP_iceil >> 2;\
        HILLOOP_iceil |= HILLOOP_iceil >> 4;\
        HILLOOP_iceil |= HILLOOP_iceil >> 8;\
        HILLOOP_iceil |= HILLOOP_iceil >> 16;\
        HILLOOP_iceil++;\
        if (HILLOOP_iceil < 8)\
            HILLOOP_iceil = 8;\
        int HILLOOP_jceil = HILLOOP_jdiff - 1;\
        HILLOOP_jceil |= HILLOOP_jceil >> 1;\
        HILLOOP_jceil |= HILLOOP_jceil >> 2;\
        HILLOOP_jceil |= HILLOOP_jceil >> 4;\
        HILLOOP_jceil |= HILLOOP_jceil >> 8;\
        HILLOOP_jceil |= HILLOOP_jceil >> 16;\
        HILLOOP_jceil++;\
        if (HILLOOP_jceil < 8)\
            HILLOOP_jceil = 8;\
        while ((i) < HILLOOP_imax && (j) < HILLOOP_jmax) {\
            int HILLOOP_icur, HILLOOP_jcur, HILLOOP_c, HILLOOP_base, HILLOOP_icase, HILLOOP_jcase;\
            unsigned long long HILLOOP_stop, HILLOOP_hilbert = 0ull;\
            if (HILLOOP_idiff > HILLOOP_jdiff) {\
                HILLOOP_jcur = HILLOOP_jdiff;\
                HILLOOP_icur = HILLOOP_imax - (i);\
                HILLOOP_base = HILLOOP_jceil / 8;\
                HILLOOP_stop = (unsigned long long) HILLOOP_base * (unsigned long long) HILLOOP_base;\
                if (HILLOOP_icur >= 2 * HILLOOP_jceil) {\
                    HILLOOP_icur = HILLOOP_jceil;\
                    HILLOOP_c = 3;\
                    HILLOOP_icase = 0;\
                    HILLOOP_jcase = HILLOOP_jcur & 1;\
                } else if (HILLOOP_icur > HILLOOP_jceil) {\
                    HILLOOP_icur = (HILLOOP_icur + 3) / 4 * 2;\
                    HILLOOP_c = 3;\
                    HILLOOP_icase = 0;\
                    HILLOOP_jcase = HILLOOP_jcur & 1;\
                } else {\
                    HILLOOP_jcase = HILLOOP_jcur & 1;\
                    HILLOOP_icase = HILLOOP_icur & 1;\
                    HILLOOP_c = 3 - (HILLOOP_icase & ~HILLOOP_jcase);\
                    HILLOOP_icase += HILLOOP_icase & HILLOOP_jcase;\
                }\
            } else {\
                HILLOOP_icur = HILLOOP_idiff;\
                HILLOOP_jcur = HILLOOP_jmax - (j);\
                HILLOOP_base = HILLOOP_iceil / 8;\
                HILLOOP_stop = (unsigned long long)HILLOOP_base * (unsigned long long)HILLOOP_base;\
                if (HILLOOP_jcur >= 2 * HILLOOP_iceil) {\
                    HILLOOP_jcur = HILLOOP_iceil;\
                    HILLOOP_c = 2;\
                    HILLOOP_jcase = 0;\
                    HILLOOP_icase = HILLOOP_icur & 1;\
                } else if (HILLOOP_jcur > HILLOOP_iceil) {\
                    HILLOOP_jcur = (HILLOOP_jcur + 3) / 4 * 2;\
                    HILLOOP_c = 2;\
                    HILLOOP_jcase = 0;\
                    HILLOOP_icase = HILLOOP_icur & 1;\
                } else {\
                    HILLOOP_jcase = HILLOOP_jcur & 1;\
                    HILLOOP_icase = HILLOOP_icur & 1;\
                    HILLOOP_c = 2 + (HILLOOP_jcase & ~HILLOOP_icase);\
                    HILLOOP_jcase += HILLOOP_icase & HILLOOP_jcase;\
                }\
            }\
            int HILLOOP_i57 = 5 + 2 * (HILLOOP_icur > HILLOOP_base * 6);\
            if(HILLOOP_icur<6) HILLOOP_i57 = HILLOOP_icur ;\
            int HILLOOP_j57 = 5 + 2 * (HILLOOP_jcur > HILLOOP_base * 6);\
            if(HILLOOP_jcur<6) HILLOOP_j57 = HILLOOP_jcur ;\
            int HILLOOP_isize, HILLOOP_jsize;\
            HILLOOP_c ^= ((HILLOOP_base & 0x55555555) == 0) ;\
            int HILLOOP_I = 0;\
            int HILLOOP_J = 0;\
            while (HILLOOP_hilbert < HILLOOP_stop) {\
                if (HILLOOP_icase)\
                    if (HILLOOP_icase == 2) {\
                        int HILLOOP_v = HILLOOP_base - 1 - HILLOOP_J;\
                        HILLOOP_v |= HILLOOP_v >> 1;\
                        HILLOOP_v |= HILLOOP_v >> 2;\
                        HILLOOP_v |= HILLOOP_v >> 4;\
                        HILLOOP_v |= HILLOOP_v >> 8;\
                        HILLOOP_v |= HILLOOP_v >> 16;\
                        HILLOOP_v = (HILLOOP_v + 1) / 2;\
                        HILLOOP_isize = (HILLOOP_I == HILLOOP_v ? HILLOOP_i57 : (HILLOOP_I < HILLOOP_v ?\
                                (unsigned long long)(HILLOOP_I + 1) * (unsigned long long)(HILLOOP_icur-HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2\
                                - (unsigned long long)HILLOOP_I * (unsigned long long)(HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2 :\
                                (unsigned long long)HILLOOP_I * (unsigned long long)(HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2\
                                - (unsigned long long)(HILLOOP_I - 1) * (unsigned long long)(HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2));\
                    } else HILLOOP_isize = HILLOOP_I == HILLOOP_base - 1 ? HILLOOP_i57 :\
                            (unsigned long long)(HILLOOP_I + 1) * (unsigned long long)(HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2\
                        - (unsigned long long)HILLOOP_I * (unsigned long long)(HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2;\
                else HILLOOP_isize = (unsigned long long)(HILLOOP_I + 1) * (unsigned long long)HILLOOP_icur / HILLOOP_base / 2 * 2\
                        - (unsigned long long)HILLOOP_I * (unsigned long long)HILLOOP_icur / HILLOOP_base / 2 * 2;\
                if (HILLOOP_jcase)\
                    if (HILLOOP_jcase == 2) {\
                        int HILLOOP_v = HILLOOP_base - 1 - HILLOOP_I;\
                        HILLOOP_v |= HILLOOP_v >> 1;\
                        HILLOOP_v |= HILLOOP_v >> 2;\
                        HILLOOP_v |= HILLOOP_v >> 4;\
                        HILLOOP_v |= HILLOOP_v >> 8;\
                        HILLOOP_v |= HILLOOP_v >> 16;\
                        HILLOOP_v = (HILLOOP_v + 1) / 2;\
                        HILLOOP_jsize = (HILLOOP_J == HILLOOP_v ? HILLOOP_j57 : (HILLOOP_J < HILLOOP_v ?\
                                (unsigned long long)(HILLOOP_J + 1) * (unsigned long long)(HILLOOP_jcur-HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2\
                                - (unsigned long long)HILLOOP_J * (unsigned long long)(HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2 :\
                                (unsigned long long)HILLOOP_J * (unsigned long long)(HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2\
                                - (unsigned long long)(HILLOOP_J - 1) * (unsigned long long)(HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2));\
                    } else HILLOOP_jsize = HILLOOP_J == HILLOOP_base - 1 ? HILLOOP_j57 :\
                            (unsigned long long)(HILLOOP_J + 1) * (unsigned long long)(HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2\
                        - (unsigned long long)HILLOOP_J * (unsigned long long)(HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2;\
                else HILLOOP_jsize = (unsigned long long)(HILLOOP_J + 1) * (unsigned long long)HILLOOP_jcur / HILLOOP_base / 2 * 2\
                        - (unsigned long long)HILLOOP_J * (unsigned long long)HILLOOP_jcur / HILLOOP_base / 2 * 2;\
                if(HILLOOP_isize>8||HILLOOP_jsize>8||HILLOOP_isize<0 || HILLOOP_jsize<0)exit(-1);\
                unsigned long long HILLOOP_nanoH = HILLOOP_nanoprog[HILLOOP_isize][HILLOOP_jsize][HILLOOP_c][0] ;\
                unsigned long long HILLOOP_nanoL = HILLOOP_nanoprog[HILLOOP_isize][HILLOOP_jsize][HILLOOP_c][1] ;\
                for(;;){


#define FUR_HILBERT_END(i,j)\
                    if(HILLOOP_nanoL == 1)\
                        break;\
                    int HILLOOP_d = HILLOOP_nanoL & 1 | HILLOOP_nanoH & 2 ;\
                    (i) += (HILLOOP_d - 2) % 2 ;\
                    (j) += (HILLOOP_d - 1) % 2 ;\
                    HILLOOP_nanoL /= 2 ;\
                    HILLOOP_nanoH /= 2 ;\
                }\
                HILLOOP_hilbert++;\
                unsigned long long HILLOOP_l = HILLOOP_hilbert & -HILLOOP_hilbert;\
                HILLOOP_l = (HILLOOP_l + HILLOOP_l/2) & 0x5555555555555555ull;\
                int HILLOOP_a = (HILLOOP_hilbert / HILLOOP_l) & 3 ;\
                int HILLOOP_isOdd = (HILLOOP_l & 0xCCCCCCCCCCCCCCCCull) != 0;\
                HILLOOP_c ^= 3 * ((HILLOOP_a==3) != HILLOOP_isOdd);\
                HILLOOP_I += (HILLOOP_c - 2) % 2 ;\
                (i) += (HILLOOP_c - 2) % 2 ;\
                HILLOOP_J += (HILLOOP_c - 1) % 2 ;\
                (j) += (HILLOOP_c - 1) % 2 ;\
                HILLOOP_c ^= ((HILLOOP_a==1) != HILLOOP_isOdd);\
            }\
        }\
    }\
}

/**************
unsigned long long HILLOOP_maxhil = (HILLOOP_imax - (i) | HILLOOP_jmax - (j)) >> 3;
HILLOOP_maxhil |= HILLOOP_maxhil >> 1; HILLOOP_maxhil |= HILLOOP_maxhil >> 2;
HILLOOP_maxhil |= HILLOOP_maxhil >> 4; HILLOOP_maxhil |= HILLOOP_maxhil >> 8;
HILLOOP_maxhil = (HILLOOP_maxhil | (HILLOOP_maxhil >> 16)) + 1;
HILLOOP_maxhil *= HILLOOP_maxhil;
for (unsigned long long HILLOOP_hilbert = 0ull ; HILLOOP_hilbert < HILLOOP_maxhil ; ){
    // Determine cell size to jump over
    // if so, do jump over
    // otherwise, do micro-program
    // determine level of cell were leaving/entering
    // make the corresponding movement

}

**************/

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

#define FGF_ALT_FOR(i, j, imax, jmax, cond) {\
    register union {\
        double HILLOOP_hflo;\
        unsigned long long HILLOOP_hint;\
    };\
    int HILLOOP_imax = (imax);\
    int HILLOOP_jmax = (jmax);\
    unsigned long long HILLOOP_hilbert = 0ull;\
    int HILLOOP_c = 2;\
    int HILLOOP_nextlevel ;\
    int HILLOOP_maxext = imax - (i) > jmax - (j) ? imax - (i) : jmax-(j) ;\
    for (HILLOOP_nextlevel = 0; HILLOOP_nextlevel < 32; HILLOOP_nextlevel++)\
        if (HILLOOP_maxext <= 1 << HILLOOP_nextlevel)\
            break;\
    int HILLOOP_maxvals = 1<<HILLOOP_nextlevel ;\
    while ((i)<HILLOOP_maxvals && (j) < HILLOOP_maxvals) {\
        if ((i)<HILLOOP_imax && (j)<HILLOOP_jmax && (cond))

#define FGF_ALT_END(i, j, cond)\
        HILLOOP_hilbert++;\
        HILLOOP_hflo = HILLOOP_hilbert & (-HILLOOP_hilbert);\
        int HILLOOP_level = (HILLOOP_hint - 4607182418800017408LL) >> 53;\
        int FURHIL_clevel ;\
        int FURHIL_lb0 = (i) ;\
        int FURHIL_ub0 = (i) + 1 ;\
        int FURHIL_lb1 = (j) ;\
        int FURHIL_ub1 = (j) + 1 ;\
        if (HILLOOP_c & 2)\
            for(FURHIL_clevel = 0 ; FURHIL_clevel < HILLOOP_nextlevel ; FURHIL_clevel++){\
                FURHIL_ub0 += (1<<FURHIL_clevel) ;\
                FURHIL_ub1 += (1<<FURHIL_clevel) ;\
                if (FURHIL_lb0 < HILLOOP_imax && FURHIL_lb1 < HILLOOP_jmax && (cond))\
                    break ;\
            }\
        else\
            for(FURHIL_clevel = 0 ; FURHIL_clevel < HILLOOP_nextlevel ; FURHIL_clevel++){\
                FURHIL_lb0 -= (1<<FURHIL_clevel) ;\
                FURHIL_lb1 -= (1<<FURHIL_clevel) ;\
                if (FURHIL_lb0 < HILLOOP_imax && FURHIL_lb1 < HILLOOP_jmax && (cond))\
                    break ;\
            }\
        if (FURHIL_clevel != 0) {\
            HILLOOP_c ^= FURHIL_clevel & 1;\
            (j) += ((HILLOOP_c - 1) % 2) * ((1 << FURHIL_clevel) - 1);\
            (i) += ((HILLOOP_c - 2) % 2) * ((1 << FURHIL_clevel) - 1);\
            HILLOOP_hilbert += (1 << (FURHIL_clevel * 2)) - 2;\
            HILLOOP_c ^= 3*(FURHIL_clevel & 1);\
        } else {\
            int HILLOOP_test = HILLOOP_level & 1;\
            int HILLOOP_action = (HILLOOP_hilbert >> (2 * HILLOOP_level)) & 3;\
            HILLOOP_c ^= 3 * (HILLOOP_test != (HILLOOP_action == 3));\
            (j) += (HILLOOP_c - 1) % 2;\
            (i) += (HILLOOP_c - 2) % 2;\
            HILLOOP_c ^= (HILLOOP_test != (HILLOOP_action == 1));\
        }\
        HILLOOP_nextlevel = HILLOOP_level ;\
    }\
}

//            printf("\t\tJUMP %d + %d =: %d\n", HILLOOP_hilbert, (1 << (FURHIL_clevel * 2)) - 2, HILLOOP_hilbert + (1 << (FURHIL_clevel * 2)) - 2);\

#endif	/* HILLOOP_H */
