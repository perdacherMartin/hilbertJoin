
#include "hilloop.h"

const long long HILLOOP_programList[][4] = {
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

const unsigned long long HILLOOP_mortonListI[] = {
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

const unsigned long long HILLOOP_mortonListJ[] = {
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


const unsigned long long HILLOOP_asuzListI[][4] = {
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

const unsigned long long HILLOOP_asuzListJ[][4] = {
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


const unsigned long long HILLOOP_nanoprog[9][9][4][2] = { {
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

const unsigned long long HILLOOP_fgfnano[4][2] = {
    {0x8ef131645ba46431ull, 0x710ece9bcef13164ull},
    {0x9ba464310ef13164ull, 0x645b9bce9ba46431ull},
    {0xa45b9bcef10ece9bull, 0x5ba46431645b9bceull},
    {0xb10ece9ba45b9bceull, 0x4ef13164310ece9bull}
} ;
