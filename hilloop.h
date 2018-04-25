/* 
 * File:   hilloop.h
 * Author: boehm
 *
 * Created on 11. März 2016, 10:49
 */

#ifndef HILLOOP_H
#define	HILLOOP_H

long long HILLOOP_programList[][4] = {
    {0LL, 0LL, 0LL, 0LL}, // 0 x 0                                                      [0]
    {0LL, 0LL, 0LL, 0LL}, // 0 x 1                                                      [1]
    {0LL, 0LL, 0LL, 0LL}, // 0 x 2                                                      [2]
    {0LL, 0LL, 0LL, 0LL}, // 0 x 3                                                      [3]
    {0LL, 0LL, 0LL, 0LL}, // 0 x 4                                                      [4]
    {0LL, 0LL, 0LL, 0LL}, // 1 x 0                                                      [5]
    {0LL, 0LL, 0LL, 0LL}, // 1 x 1                                                      [6]
    {0xeLL, 0xeLL, 0xeLL, 0xeLL}, // 1 x 2                                              [7]
    {0xeeLL, 0xeeLL, 0xeeLL, 0xeeLL}, // 1 x 3                                          [8]
    {0xeeeLL, 0xeeeLL, 0xeeeLL, 0xeeeLL}, // 1 x 4                                      [9]
    {0LL, 0LL, 0LL, 0LL}, // 2 x 0                                                      [10]
    {0xbLL, 0xbLL, 0xbLL, 0xbLL}, // 2 x 1                                              [11]
    {0xb69LL, 0xe96LL, 0x9ebLL, 0x6beLL}, // i x j = 2 x 2 ;                            [12]
    {0xb5b69LL, 0xee966LL, 0x9f9ebLL, 0x66beeLL}, // i x j = 2 x 3 ;                    [13]
    {0xb696b69LL, 0xeee9666LL, 0x9ebe9ebLL, 0x666beeeLL}, // i x j = 2 x 4 ;            [14]
    {0LL, 0LL, 0LL, 0LL}, // 3 x 0                                                      [15]
    {0xbbLL, 0xbbLL, 0xbbLL, 0xbbLL}, // 3 x 1                                          [16]
    {0xbb699LL, 0xe5e96LL, 0x99ebbLL, 0x6f6beLL}, // i x j = 3 x 2 ;                    [17]
    {0xbb669e96LL, 0xee996b69LL, 0x99eeb6beLL, 0x66bbe9ebLL}, // i x j = 3 x 3 ;        [18]
    {0xbb6996bb699LL, 0x9f9ebe99666LL, 0x99ebbe99ebbLL, 0xb5b696bbeeeLL}, // i x j = 3 x 4 ; [19]
    {0LL, 0LL, 0LL, 0LL}, // 4 x 0                                                      [20]
    {0xbbbLL, 0xbbbLL, 0xbbbLL, 0xbbbLL}, // 4 x 1                                      [21]
    {0xbbb6999LL, 0xe969e96LL, 0x999ebbbLL, 0x6beb6beLL}, // i x j = 4 x 2 ;            [22]
    {0x6f6beb66999LL, 0xee9669ee966LL, 0xe5e969eebbbLL, 0x66beeb66beeLL}, // i x j = 4 x 3 ; [23]
    {0x6bebb696b699e96LL, 0x9ebee969e966b69LL, 0xe9699ebe9ebb6beLL, 0xb6966beb6bee9ebLL} // i x j = 4 x 4 ; [24]
};

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

# define HILLOOP_END(i,j)\
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

unsigned long long HILLOOP_mortonListI[] = {
    0ull,               // 0 x 0: 0
    0ull,               // 0 x 1: 1
    0ull,               // 0 x 2: 2
    0ull,               // 0 x 3: 3
    0ull,               // 0 x 4: 4
    0ull,               // 1 x 0: 5
    05ull,              // 1 x 1: 6
    055ull,             // 1 x 2: 7
    0555ull,            // 1 x 3: 8
    05555ull,           // 1 x 4: 9
    0ull,               // 2 x 0: 10
    056ull,             // 2 x 1: 11
    05565ull,           // 2 x 2: 12
    0564565ull,         // 2 x 3: 13
    055654565ull,       // 2 x 4: 14
    0ull,               // 3 x 0: 15
    0566ull,            // 3 x 1: 16
    0556565ull,         // 3 x 2: 17
    0555664565ull,      // 3 x 3: 18
    0556565356565ull,   // 3 x 4: 19
    0ull,               // 4 x 0: 20
    05666ull,           // 4 x 1: 21
    055656565ull,       // 4 x 2: 22
    0564565664565ull,   // 4 x 3: 23
    05565456565654565ull// 4 x 4: 24
};

unsigned long long HILLOOP_mortonListJ[] = {
    0ull,               // 0 x 0: 0
    0ull,               // 0 x 1: 1
    0ull,               // 0 x 2: 2
    0ull,               // 0 x 3: 3
    0ull,               // 0 x 4: 4
    0ull,               // 1 x 0: 5
    05ull,              // 1 x 1: 6
    056ull,             // 1 x 2: 7
    0566ull,            // 1 x 3: 8
    05666ull,           // 1 x 4: 9
    0ull,               // 2 x 0: 10
    055ull,             // 2 x 1: 11
    05646ull,           // 2 x 2: 12
    0556646ull,         // 2 x 3: 13
    056466646ull,       // 2 x 4: 14
    0ull,               // 3 x 0: 15
    0555ull,            // 3 x 1: 16
    0564646ull,         // 3 x 2: 17
    0566356646ull,      // 3 x 3: 18
    0564646664646ull,   // 3 x 4: 19
    0ull,               // 3 x 0: 20
    05555ull,           // 4 x 1: 21
    056464646ull,       // 4 x 2: 22
    0556646356646ull,   // 4 x 3: 23
    05646664626466646ull// 4 x 4: 24
};


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

# define ZORDLOOP_END(i,j)\
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

unsigned long long HILLOOP_asuzListI[][4] = {
    {0ull,0ull,0ull,0ull},                              // 0 x 0: 0
    {0ull,0ull,0ull,0ull},                              // 0 x 1: 1
    {0ull,0ull,0ull,0ull},                              // 0 x 2: 2
    {0ull,0ull,0ull,0ull},                              // 0 x 3: 3
    {0ull,0ull,0ull,0ull},                              // 0 x 4: 4
    {0ull,0ull,0ull,0ull},                              // 1 x 0: 5
    {05ull,05ull,05ull,05ull},                          // 1 x 1: 6
    {055ull,055ull,055ull,055ull},                      // 1 x 2: 7
    {0555ull,0555ull,0555ull,0555ull},                  // 1 x 3: 8
    {05555ull,05555ull,05555ull,05555ull},                                                // 1 x 4: 9
    {0ull,0ull,0ull,0ull},                                                                // 2 x 0: 10
    {056ull,056ull,056ull,056ull},                                                        // 2 x 1: 11
    {05654ull,05545ull,05456ull,05565ull},                                                // 2 x 2: 12
    {0564654ull,0555455ull,0546456ull,0555655ull},                                        // 2 x 3: 13
    {056545654ull,055456545ull,054565456ull,055654565ull},                                // 2 x 4: 14
    {0ull,0ull,0ull,0ull},                                                                // 3 x 0: 15
    {0566ull,0566ull,0566ull,0566ull},                                                    // 3 x 1: 16
    {0566544ull,0554545ull,0544566ull,0556565ull},                                        // 3 x 2: 17
    {0556565544ull,0555455455ull,0554545566ull,0555655655ull},                            // 3 x 3: 18
    {0566544566544ull,0544566554545ull,0544566544566ull,0566544556565ull},                // 3 x 4: 19
    {0ull,0ull,0ull,0ull},                                                                // 4 x 0: 20
    {05666ull,05666ull,05666ull,05666ull},                                                // 4 x 1: 21
    {056665444ull,055454545ull,054445666ull,055656565ull},                                // 4 x 2: 22
    {0556565655444ull,0555455455455ull,0554545455666ull,0555655655655ull},                // 4 x 3: 23
    {05565665456544545ull,05456554545455654ull,05545445654566565ull,05654556565655456ull} // 4 x 4: 24
};

unsigned long long HILLOOP_asuzListJ[][4] = {
    {0ull,0ull,0ull,0ull},                              // 0 x 0: 0
    {0ull,0ull,0ull,0ull},                              // 0 x 1: 1
    {0ull,0ull,0ull,0ull},                              // 0 x 2: 2
    {0ull,0ull,0ull,0ull},                              // 0 x 3: 3
    {0ull,0ull,0ull,0ull},                              // 0 x 4: 4
    {0ull,0ull,0ull,0ull},                              // 1 x 0: 5
    {05ull,05ull,05ull,05ull},                          // 1 x 1: 6
    {056ull,056ull,056ull,056ull},                      // 1 x 2: 7
    {0566ull,0566ull,0566ull,0566ull},                  // 1 x 3: 8
    {05666ull,05666ull,05666ull,05666ull},              // 1 x 4: 9
    {0ull,0ull,0ull,0ull},                              // 2 x 0: 10
    {055ull,055ull,055ull,055ull},                      // 2 x 1: 11
    {05565ull,05646ull,05565ull,05646ull},              // 2 x 2: 12
    {0556565ull,0566366ull,0556565ull,0566366ull},      // 2 x 3: 13
    {055656565ull,056466646ull,055656565ull,056466646ull},                              // 2 x 4: 14
    {0ull,0ull,0ull,0ull},                                                              // 3 x 0: 15
    {0555ull,0555ull,0555ull,0555ull},                                                  // 3 x 1: 16
    {0555655ull,0564646ull,0555655ull,0564646ull},                                      // 3 x 2: 17
    {0564646655ull,0566366366ull,0564646655ull,0566366366ull},                          // 3 x 3: 18
    {0555655655655ull,0555655664646ull,0555655655655ull,0555655664646ull},              // 3 x 4: 19
    {0ull,0ull,0ull,0ull},                                                              // 4 x 0: 20
    {05555ull,05555ull,05555ull,05555ull},                                              // 4 x 1: 21
    {055556555ull,056464646ull,055556555ull,056464646ull},                              // 4 x 2: 22
    {0564646466555ull,0566366366366ull,0564646466555ull,0566366366366ull},              // 4 x 3: 23
    {05646456565654646ull,05565664626466565ull,05646456565654646ull,05565664626466565ull}// 4 x 4: 24
};

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

# define ASUZLOOP_END(i,j)\
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

# define BIGHLOOP_END(i,j)\
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

# define PEANOLOOP_END(i,j)\
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

# define BIGZLOOP_END(i,j)\
        }\
    }\
}\

unsigned long long HILLOOP_nanoprog[9][9][4][2] = { {
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 0
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 1
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 2
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 3
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 4
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 5
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 6
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 0 x 7
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} } }, // 0 x 8
  { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 1 x 0
    { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 1 x 1
    { {0x0ull, 0x2ull}, {0x0ull, 0x2ull}, {0x2ull, 0x2ull}, {0x2ull, 0x2ull} }, // 1 x 2
    { {0x0ull, 0x4ull}, {0x0ull, 0x4ull}, {0x6ull, 0x4ull}, {0x6ull, 0x4ull} }, // 1 x 3
    { {0x0ull, 0x8ull}, {0x0ull, 0x8ull}, {0xeull, 0x8ull}, {0xeull, 0x8ull} }, // 1 x 4
    { {0x0ull, 0x10ull}, {0x0ull, 0x10ull}, {0x1eull, 0x10ull}, {0x1eull, 0x10ull} }, // 1 x 5
    { {0x0ull, 0x20ull}, {0x0ull, 0x20ull}, {0x3eull, 0x20ull}, {0x3eull, 0x20ull} }, // 1 x 6
    { {0x0ull, 0x40ull}, {0x0ull, 0x40ull}, {0x7eull, 0x40ull}, {0x7eull, 0x40ull} }, // 1 x 7
    { {0x0ull, 0x80ull}, {0x0ull, 0x80ull}, {0xfeull, 0x80ull}, {0xfeull, 0x80ull} } }, // 1 x 8
  { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 2 x 0
    { {0x0ull, 0x3ull}, {0x0ull, 0x3ull}, {0x2ull, 0x3ull}, {0x2ull, 0x3ull} }, // 2 x 1
    { {0x8ull, 0xdull}, {0x8ull, 0xaull}, {0x6ull, 0xdull}, {0x6ull, 0xaull} }, // 2 x 2
    { {0x0ull, 0x1ull}, {0x30ull, 0x24ull}, {0x0ull, 0x1ull}, {0xeull, 0x24ull} }, // 2 x 3
    { {0x88ull, 0xd5ull}, {0xe0ull, 0x88ull}, {0x76ull, 0xd5ull}, {0x1eull, 0x88ull} }, // 2 x 4
    { {0x0ull, 0x1ull}, {0x3c0ull, 0x210ull}, {0x0ull, 0x1ull}, {0x3eull, 0x210ull} }, // 2 x 5
    { {0x888ull, 0xd55ull}, {0xf80ull, 0x820ull}, {0x776ull, 0xd55ull}, {0x7eull, 0x820ull} }, // 2 x 6
    { {0x0ull, 0x1ull}, {0x3f00ull, 0x2040ull}, {0x0ull, 0x1ull}, {0xfeull, 0x2040ull} }, // 2 x 7
    { {0x8888ull, 0xd555ull}, {0xfe00ull, 0x8080ull}, {0x7776ull, 0xd555ull}, {0x1feull, 0x8080ull} } }, // 2 x 8
  { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 3 x 0
    { {0x0ull, 0x7ull}, {0x0ull, 0x7ull}, {0x6ull, 0x7ull}, {0x6ull, 0x7ull} }, // 3 x 1
    { {0x30ull, 0x3bull}, {0x0ull, 0x1ull}, {0xeull, 0x3bull}, {0x0ull, 0x1ull} }, // 3 x 2
    { {0x188ull, 0x1caull}, {0x188ull, 0x135ull}, {0x76ull, 0x1caull}, {0x76ull, 0x135ull} }, // 3 x 3
    { {0x708ull, 0xa8aull}, {0x0ull, 0x1ull}, {0x8f6ull, 0xa8aull}, {0x0ull, 0x1ull} }, // 3 x 4
    { {0x6230ull, 0x729bull}, {0x3e20ull, 0x68d4ull}, {0x1dceull, 0x729bull}, {0x41deull, 0x68d4ull} }, // 3 x 5
    { {0x30c30ull, 0x3b6dbull}, {0x0ull, 0x1ull}, {0xf3ceull, 0x3b6dbull}, {0x0ull, 0x1ull} }, // 3 x 6
    { {0x188c30ull, 0x1ca6dbull}, {0x1be208ull, 0x128d45ull}, {0x773ceull, 0x1ca6dbull}, {0x41df6ull, 0x128d45ull} }, // 3 x 7
    { {0xc30c30ull, 0xedb6dbull}, {0x0ull, 0x1ull}, {0x3cf3ceull, 0xedb6dbull}, {0x0ull, 0x1ull} } }, // 3 x 8
  { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 4 x 0
    { {0x0ull, 0xfull}, {0x0ull, 0xfull}, {0xeull, 0xfull}, {0xeull, 0xfull} }, // 4 x 1
    { {0xe0ull, 0xf7ull}, {0x88ull, 0xaaull}, {0x1eull, 0xf7ull}, {0x76ull, 0xaaull} }, // 4 x 2
    { {0x0ull, 0x1ull}, {0x708ull, 0xd75ull}, {0x0ull, 0x1ull}, {0x8f6ull, 0xd75ull} }, // 4 x 3
    { {0x7888ull, 0xad5aull}, {0x7888ull, 0xd2a5ull}, {0x8776ull, 0xad5aull}, {0x8776ull, 0xd2a5ull} }, // 4 x 4
    { {0x0ull, 0x1ull}, {0x7c308ull, 0xd1245ull}, {0x0ull, 0x1ull}, {0x83cf6ull, 0xd1245ull} }, // 4 x 5
    { {0x1f8888ull, 0x88d55aull}, {0x778888ull, 0xd52a55ull}, {0xe07776ull, 0x88d55aull}, {0x887776ull, 0xd52a55ull} }, // 4 x 6
    { {0x0ull, 0x1ull}, {0x77c3088ull, 0xd512455ull}, {0x0ull, 0x1ull}, {0x883cf76ull, 0xd512455ull} }, // 4 x 7
    { {0x78887888ull, 0xad5a2d5aull}, {0x77788888ull, 0xd552a555ull}, {0x87778776ull, 0xad5a2d5aull}, {0x88877776ull, 0xd552a555ull} } }, // 4 x 8
  { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 5 x 0
    { {0x0ull, 0x1full}, {0x0ull, 0x1full}, {0x1eull, 0x1full}, {0x1eull, 0x1full} }, // 5 x 1
    { {0x3c0ull, 0x3efull}, {0x0ull, 0x1ull}, {0x3eull, 0x3efull}, {0x0ull, 0x1ull} }, // 5 x 2
    { {0x3e20ull, 0x572bull}, {0x6230ull, 0x4d64ull}, {0x41deull, 0x572bull}, {0x1dceull, 0x4d64ull} }, // 5 x 3
    { {0x7c308ull, 0xaedbaull}, {0x0ull, 0x1ull}, {0x83cf6ull, 0xaedbaull}, {0x0ull, 0x1ull} }, // 5 x 4
    { {0x7e2308ull, 0x12729baull}, {0x7e2308ull, 0x1d8d645ull}, {0x181dcf6ull, 0x12729baull}, {0x181dcf6ull, 0x1d8d645ull} }, // 5 x 5
    { {0x7f0c308ull, 0x223b6dbaull}, {0x0ull, 0x1ull}, {0x380f3cf6ull, 0x223b6dbaull}, {0x0ull, 0x1ull} }, // 5 x 6
    { {0x1f88c30e0ull, 0x49ca6db88ull}, {0x1e7e23088ull, 0x76d8d6455ull}, {0x60773cf1eull, 0x49ca6db88ull}, {0x6181dcf76ull, 0x76d8d6455ull} }, // 5 x 7
    { {0x7c3087c308ull, 0xaedba2edbaull}, {0x0ull, 0x1ull}, {0x83cf783cf6ull, 0xaedba2edbaull}, {0x0ull, 0x1ull} } }, // 5 x 8
  { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 6 x 0
    { {0x0ull, 0x3full}, {0x0ull, 0x3full}, {0x3eull, 0x3full}, {0x3eull, 0x3full} }, // 6 x 1
    { {0xf80ull, 0xfdfull}, {0x888ull, 0xaaaull}, {0x7eull, 0xfdfull}, {0x776ull, 0xaaaull} }, // 6 x 2
    { {0x0ull, 0x1ull}, {0x30c30ull, 0x24924ull}, {0x0ull, 0x1ull}, {0xf3ceull, 0x24924ull} }, // 6 x 3
    { {0x778888ull, 0xaad5aaull}, {0x1f8888ull, 0xf72aa5ull}, {0x887776ull, 0xaad5aaull}, {0xe07776ull, 0xf72aa5ull} }, // 6 x 4
    { {0x0ull, 0x1ull}, {0x7f0c308ull, 0x3dc49245ull}, {0x0ull, 0x1ull}, {0x380f3cf6ull, 0x3dc49245ull} }, // 6 x 5
    { {0x1f7888e08ull, 0x88ad5a77aull}, {0x1f7888e08ull, 0xf752a5885ull}, {0xe087771f6ull, 0x88ad5a77aull}, {0xe087771f6ull, 0xf752a5885ull} }, // 6 x 6
    { {0x0ull, 0x1ull}, {0x21ddf0c3088ull, 0x2b568492455ull}, {0x0ull, 0x1ull}, {0x1e220f3cf76ull, 0x2b568492455ull} }, // 6 x 7
    { {0x877788887888ull, 0xd2a5d555d2a5ull}, {0x87777888e088ull, 0xad5a52a58855ull}, {0x788877778776ull, 0xd2a5d555d2a5ull}, {0x788887771f76ull, 0xad5a52a58855ull} } }, // 6 x 8
  { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 7 x 0
    { {0x0ull, 0x7full}, {0x0ull, 0x7full}, {0x7eull, 0x7full}, {0x7eull, 0x7full} }, // 7 x 1
    { {0x3f00ull, 0x3fbfull}, {0x0ull, 0x1ull}, {0xfeull, 0x3fbfull}, {0x0ull, 0x1ull} }, // 7 x 2
    { {0x1be208ull, 0x1d72baull}, {0x188c30ull, 0x135924ull}, {0x41df6ull, 0x1d72baull}, {0x773ceull, 0x135924ull} }, // 7 x 3
    { {0x77c3088ull, 0xaaedbaaull}, {0x0ull, 0x1ull}, {0x883cf76ull, 0xaaedbaaull}, {0x0ull, 0x1ull} }, // 7 x 4
    { {0x1e7e23088ull, 0x492729baaull}, {0x1f88c30e0ull, 0x763592477ull}, {0x6181dcf76ull, 0x492729baaull}, {0x60773cf1eull, 0x763592477ull} }, // 7 x 5
    { {0x21ddf0c3088ull, 0x34a97b6dbaaull}, {0x0ull, 0x1ull}, {0x1e220f3cf76ull, 0x34a97b6dbaaull}, {0x0ull, 0x1ull} }, // 7 x 6
    { {0x79f88c307888ull, 0x1249ca6dbd2a5ull}, {0x79f88c307888ull, 0x1db6359242d5aull}, {0x1860773cf8776ull, 0x1249ca6dbd2a5ull}, {0x1860773cf8776ull, 0x1db6359242d5aull} }, // 7 x 7
    { {0x8777c30c307888ull, 0xd2a5edb6dbd2a5ull}, {0x0ull, 0x1ull}, {0x78883cf3cf8776ull, 0xd2a5edb6dbd2a5ull}, {0x0ull, 0x1ull} } }, // 7 x 8
  { { {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull}, {0x0ull, 0x1ull} }, // 8 x 0
    { {0x0ull, 0xffull}, {0x0ull, 0xffull}, {0xfeull, 0xffull}, {0xfeull, 0xffull} }, // 8 x 1
    { {0xfe00ull, 0xff7full}, {0x8888ull, 0xaaaaull}, {0x1feull, 0xff7full}, {0x7776ull, 0xaaaaull} }, // 8 x 2
    { {0x0ull, 0x1ull}, {0xc30c30ull, 0x924924ull}, {0x0ull, 0x1ull}, {0x3cf3ceull, 0x924924ull} }, // 8 x 3
    { {0x77788888ull, 0xaaad5aaaull}, {0x78887888ull, 0xd2a5d2a5ull}, {0x88877776ull, 0xaaad5aaaull}, {0x87778776ull, 0xd2a5d2a5ull} }, // 8 x 4
    { {0x0ull, 0x1ull}, {0x7c3087c308ull, 0xd1245d1245ull}, {0x0ull, 0x1ull}, {0x83cf783cf6ull, 0xd1245d1245ull} }, // 8 x 5
    { {0x87777888e088ull, 0xd2a5ad5a77aaull}, {0x877788887888ull, 0xad5a2aaa2d5aull}, {0x788887771f76ull, 0xd2a5ad5a77aaull}, {0x788877778776ull, 0xad5a2aaa2d5aull} }, // 8 x 6
    { {0x0ull, 0x1ull}, {0x8777c30c307888ull, 0xad5a1249242d5aull}, {0x0ull, 0x1ull}, {0x78883cf3cf8776ull, 0xad5a1249242d5aull} }, // 8 x 7
    { {0x8777788878887888ull, 0xd2a5ad5a2d5ad2a5ull}, {0x8777788878887888ull, 0xad5a52a5d2a52d5aull}, {0x7888877787778776ull, 0xd2a5ad5a2d5ad2a5ull}, {0x7888877787778776ull, 0xad5a52a5d2a52d5aull} } }// 8 x 8
} ;

unsigned long long HILLOOP_fgfnano[4][2] = {
    {0x8ef131645ba46431ull, 0x710ece9bcef13164ull},
    {0x9ba464310ef13164ull, 0x645b9bce9ba46431ull},
    {0xa45b9bcef10ece9bull, 0x5ba46431645b9bceull},
    {0xb10ece9ba45b9bceull, 0x4ef13164310ece9bull}
} ;

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
    
                    
# define FUR_HILBERT_END(i,j)\
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

# define FGF_HILBERT_FOR(i, j, imax, jmax, cond, cond2) {\
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

# define FGF_HILBERT_END(i, j)\
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

# define FGFN_HILBERT_FOR(i, j, imax, jmax, cond, cond2) {\
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

# define FGFN_HILBERT_END(i, j)\
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

# define FGFC_HILBERT_FOR(i, j, imax, jmax, cond, cond2, cond3) {\
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

# define FGFC_HILBERT_ELSE(i, j)\
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

# define FGFC_HILBERT_END(i, j)\
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

# define FGF_ALT_FOR(i, j, imax, jmax, cond) {\
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

# define FGF_ALT_END(i, j, cond)\
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

