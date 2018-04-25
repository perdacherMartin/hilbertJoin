/* update version:
 * Email vom 5. Oktber 2017
 */

#include <stdio.h>
#include <string.h>

#include "measure/timer.h"
#include "measure/energy.h"
#include "util/dataIo.h"
#include "util/chrisutil.h"
#include "util/arguments.h"
#include "hilbertjoin/egojoin.h"


int main(int argc, char** argv) {
    size_t n = 4000000;
    size_t d = 41;
    size_t threads=64;
    double epsilon = 0.034;
    char filename[256] = "";
    bool isBinary=false;
    CUtilTimer timer;
    Hioki pmeter;
    size_t result=0l;
    int KBLOCK=8,stripes=2;

    parsing_args(argc, argv, &n, &epsilon, &d, &threads, filename, &isBinary,&KBLOCK,&stripes);

    omp_set_num_threads(threads);
    // double * array = (double*) mallocA64((n+7)/8*8 * sizeof (double) * d + 16384);
    // double * array = (double*) mallocA64(n * sizeof (double) * d + 16384);
    // double * array = (double*) ddr_alloc((n+7)/8*8 * sizeof(double) * d + 16384);
    double * array = (double*) ddr_alloc(n * sizeof (double) * d + 16384);
    // printf("alloc ok\n"); fflush(stdout);
    read_file(array, n, d, filename, isBinary);
    // printf("readfile ok\n"); fflush(stdout);
    pmeter.reset(); pmeter.start();
    timer.start();
    // test_ego_loop3(n,d,threads,epsilon,array,&result);
    // printf("start\n"); fflush(stdout);
    // test_ego_loop3_long(n,d,threads,epsilon,array,&result,stripes,KBLOCK);
    test_ego_loop3_macro(n,d,threads,epsilon,array,&result,stripes,KBLOCK);
    timer.stop();
    pmeter.stop();

    printf("%zu;%zu;%zu;%d;%d;%f;%f;%ld;%f\n", n,d,threads,epsilon,stripes,KBLOCK,timer.get_time(),result,pmeter.getWH());
    // freeA64(array);
    ddr_free(array);

    return 0;
}
