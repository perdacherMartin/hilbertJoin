
// main method for self-join
// for a two-set join see mainJoin.cpp

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#include "measure/timer.h"
#include "measure/energy.h"
#include "util/dataIo.h"
#include "util/chrisutil.h"
#include "util/arguments.h"
#include "hilbertjoin/egojoin.h"

#include <boost/lockfree/queue.hpp>
#include <boost/atomic.hpp>

#ifndef COUNT_ONLY
size_t consumer_count;

int consumer(boost::lockfree::queue<join_pair> &queue)
{
    join_pair jp;
    while (queue.pop(jp)){
        #pragma omp atomic write
            consumer_count= consumer_count +1;
        // printf("%lu-%lu\n", jp.p1, jp.p2);
    }
}
#endif

int main(int argc, char** argv) {
    size_t n = 4000000;
    size_t d = 41;
    size_t threads=64;
    double epsilon = 0.034;
    char filename[256] = "";
    bool isBinary=false;
    CUtilTimer timer, algtimer;
    Hioki pmeter;
    size_t result=0l;
    int stripes=14;
    int actdim=3;
    boost::lockfree::queue<join_pair> queue(10000);
    double sortTime=0.0, reorderTime=0.0, indexTime=0.0, watthours=0.0,totaltime=0.0,algtime=0.0;
    double loadpercent=0.0;

    parsing_args(argc, argv, &n, &epsilon, &d, filename, &isBinary, &actdim);

    stripes = ((int)pow(3,actdim) + 1) / 2;
    omp_set_num_threads(NUM_THREADS);
    int *reorder_dim=(int*) malloc ((d+8)*sizeof(int));
// #ifndef COUNT_ONLY
//     printf("COUNT_ONLY is not defined.");
// #endif
//
// #ifdef COUNT_ONLY
//     printf("COUNT_ONLY is defined.");
// #endif
    // printf("Using %d threads!\n", NUM_THREADS);

    // omp_set_num_threads(threads);
    // double * array = (double*) mallocA64((n+7)/8*8 * sizeof (double) * d + 16384);
    // double * array = (double*) mallocA64(n * sizeof (double) * d + 16384);
    // double * array = (double*) ddr_alloc((n+7)/8*8 * sizeof(double) * d + 16384);
    double * array = (double*) ddr_alloc(n * sizeof (double) * d + 16384);
    // printf("alloc ok\n"); fflush(stdout);
    read_file(array, n, d, filename, isBinary);
    // printf("readfile ok\n"); fflush(stdout);

    for ( int i=0 ; i < 10 ; i++ ){
        for ( int j=0 ; j < d ; j++ ){
            printf("%f, ", array[i*d+j]);
        }
        printf("\n");
    }

    pmeter.reset(); pmeter.start();
    timer.start();

    outputStatistics(n, d, epsilon, array, reorder_dim);
    // sampleHistograms(n, d, epsilon, array, reorder_dim);
    reorder_dimensions(n, d, array, reorder_dim);


    timer.stop();
    reorderTime = timer.get_time();
    // test_ego_loop3(n,d,threads,epsilon,array,&result);
    // printf("start\n"); fflush(stdout);
    // test_ego_loop3_long(n,d,threads,epsilon,array,&result,stripes,KBLOCK);
    algtimer.start();
#ifdef COUNT_ONLY
    test_ego_loop3_macro(n,d,epsilon,array,&result,stripes,&sortTime,&indexTime,&loadpercent);
    // test_ego_loop3_macro(n,d,threads,epsilon,array,&result,stripes,&sortTime);
#else
    // test_ego_loop3_macro_queue(n,d,threads,epsilon,array,&result,stripes,KBLOCK,queue, &sortTime);
#endif

    // test_ego_loop3_macro_queue(n,d,threads,epsilon,array,&result,stripes,KBLOCK,queue);
    // test(queue);

    algtimer.stop();
    algtime = algtimer.get_time();
    pmeter.stop();
    watthours=pmeter.getWH();

#ifndef COUNT_ONLY
    // if we materialize with a non-blocking linked list, then joincounts are zero
    #pragma omp parallel for
    for ( int i = 0 ; i < threads ; i++ ){
        consumer(queue);
    }
    // printf("overwrite result\n");
    result = consumer_count;
#endif
    double jp_per_point = (result == 0 ) ? 0 : (double)result / n ;
    // HEADER:
    // N;D;JPPP;THREADS;EPSILON;STRIPES;KBLOCK;TIME;ALGTIME;SORTTIME;INDEXTIME;REORDERTIME;COUNTS;LOADPERCENT;WH
    printf("%zu;%zu;%f;%zu;%2.14f;%d;%d;%f;%f;%f;%f;%f;%ld;%f;%f\n", n,d,jp_per_point, NUM_THREADS,epsilon,stripes,KBLOCK,algtime+reorderTime,algtime - sortTime,sortTime,indexTime,reorderTime,result,loadpercent,watthours);
    // freeA64(array);
    ddr_free(array);
    free(reorder_dim);

    return 0;
}
