
#include "papicalls.h"

PapiBin::PapiBin(){
    int num_hwcntrs;
    if ((num_hwcntrs = PAPI_num_counters()) < 0 ){
        if ((num_hwcntrs = PAPI_num_counters()) < NUM_EVENTS ){
            fprintf(stderr,"Info:: This installation does not support PAPI: %s\n", PAPI_strerror(num_hwcntrs));
            exit(1);
        }
    }

    if ((num_hwcntrs = PAPI_num_counters()) < NUM_EVENTS ){
        fprintf(stderr,"Info:: This machine does not provide sufficient hardware counters.\n");
        exit(1);
    }
}

void PapiBin::start(){
    int retval=0;

    counters.l1 = 0; // counters.l2 = 0; counters.l3 = 0;
    // counters.rtime = 0.0; counters.ptime = 0.0; counters.flops = 0; counters.mflops;

    // if ( ( retval = PAPI_flops( &counters.rtime, &counters.ptime, &counters.flops, &counters.mflops ) ) < PAPI_OK ){
    //     fprintf(stderr,"Call to PAPI_flops failed: %s\n", PAPI_strerror(retval));
    //     exit(1);
    // }
    timer.start();

    if ((retval = PAPI_start_counters(events, NUM_EVENTS)) < PAPI_OK) {
        fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(retval));
        exit(1);
    }
}

void PapiBin::stop(){
    int retval=0;
    long_long values[NUM_EVENTS];

    if ((retval = PAPI_stop_counters(values, NUM_EVENTS)) < PAPI_OK) {
        fprintf(stderr, "PAPI failed to start counters: %s\n", PAPI_strerror(retval));
        exit(1);
    }

    // if ( ( retval = PAPI_flops( &counters.rtime, &counters.ptime, &counters.flops, &counters.mflops ) ) < PAPI_OK ){
    //     fprintf(stderr,"Call to PAPI_flops failed: %s\n", PAPI_strerror(retval));
    //     exit(1);
    // }
    timer.stop();

    counters.l1 = values[0];
    counters.rtime = timer.get_time();

    // unfortunately the following do not work for xeon-phi:
    // counters.l2 = values[1];
    // counters.l3 = values[2];
}

CounterBin PapiBin::getBin(){
    return counters;
}
