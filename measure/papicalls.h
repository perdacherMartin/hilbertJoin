
#ifndef CACHE_HIERACHY_H
#define CACHE_HIERACHY_H

#include <papi.h>

#include <stdio.h>
#include <stdlib.h>

#include "timer.h"

#define NUM_EVENTS 1

struct CounterBin{
    long_long l1;
    long_long l2;
    // long_long l3;
    double rtime;
    // float ptime;
    // long_long flops;
    // float mflops;
};

class PapiBin{
public:
    PapiBin();
    void start();
    void stop();
    CounterBin getBin();
private:
    CUtilTimer timer;
    // const int NUM_EVENTS=3;
    // int events[NUM_EVENTS] = { PAPI_L1_DCM , PAPI_L2_DCM, PAPI_L3_DCM };
    int events[NUM_EVENTS] = { PAPI_L1_DCM };
    CounterBin counters;
};

#endif
