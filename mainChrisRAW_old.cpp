/*
 * File:   main.cpp
 * Author: boehm
 *
 * Created on 17. März 2017, 14:13
 */

using namespace std;
void *__gxx_personality_v0;

#ifndef NUM_THREADS
#define NUM_THREADS 64
#endif

typedef double vec __attribute__((vector_size(64), aligned(64)));
typedef double vec4 __attribute__((vector_size(32), aligned(32)));
typedef double vecu __attribute__((vector_size(64)));
typedef long long veci64 __attribute__((vector_size(64), aligned(64)));

#include <omp.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <immintrin.h>
#include <x86intrin.h>
#include <math.h>
#include "hilloop.h"

#define DBL_MAX         1.7976931348623158e+308

//#include <cstdlib>

double * mallocA64(size_t s) {
    long long adr = (long long) malloc(s + 72);
    long long adr2 = (adr + 71) & ~63;
    ((long long *) adr2)[-1] = adr;
    return (double *) adr2;
}

double * callocA64(size_t s) {
    long long adr = (long long) calloc(s + 72, 1);
    long long adr2 = (adr + 71) & ~63;
    ((long long *) adr2)[-1] = adr;
    return (double *) adr2;
}

void freeA64(void * adr) {
    free((void *) (((long long *) adr)[-1]));
}

double stopc(void) {
    static double last;
    double c = clock();
    double r = c - last;
    last = c;
    return r / CLOCKS_PER_SEC;
}

double stop(void) {
    static struct timeval last;
    struct timeval c;
    gettimeofday(&c, NULL);
    double r = (double) (c.tv_sec - last.tv_sec) + (double) (c.tv_usec - last.tv_usec) / 1000000.;
    last.tv_sec = c.tv_sec;
    last.tv_usec = c.tv_usec;
    return r;
}

typedef struct TreeEl TreeEl;
struct TreeEl {
    int n;
    int minid;
    omp_lock_t lock;
    TreeEl* left;
    TreeEl* right;
} ;
typedef struct ListEl ListEl;
struct ListEl {
    int n;
    int minid;
    omp_lock_t lock;
    ListEl* left;
} ;

TreeEl ** tentArray;
ListEl ** listArray;

void overwrite(TreeEl *smaller, TreeEl *bigger) {
    if(smaller->left == (TreeEl *)0){
        tentArray[smaller->minid] = bigger;
    } else {
        overwrite(smaller->left, bigger);
        overwrite(smaller->right, bigger);
    }
}

void integrate(TreeEl *smaller, TreeEl *bigger) {
    TreeEl *h = (TreeEl *) malloc(sizeof(TreeEl));
    h->left = bigger->left;
    h->right = bigger->right;
    h->minid = bigger->minid;
    h->n = bigger->n;
    bigger->left = h;
    bigger->right = smaller;
    bigger->n += smaller->n;
    if(smaller->minid < bigger->minid)
        bigger->minid = smaller->minid;
}

void unifyTree(int i, int j, int count) {
    if (tentArray[i] == tentArray[j])
        return;
    TreeEl * mylock1;
    TreeEl * mylock2;
    for(;;) {
        mylock1 = tentArray[i];
        omp_set_lock(&mylock1->lock);
        if (tentArray[i] == tentArray[j]) {
            omp_unset_lock(&mylock1->lock);
            return;
        }
        if (mylock1 == tentArray[i])
            break;
        omp_unset_lock(&mylock1->lock);
    }
    for (;;) {
        while (omp_test_lock(&(mylock2=tentArray[j])->lock)){
            if (mylock2 == tentArray[j])
                goto locksOkay;
            omp_unset_lock(&mylock2->lock);
        }
        do {
            omp_unset_lock(&mylock1->lock);
            mylock1 = tentArray[j];
            omp_set_lock(&mylock1->lock);
            if (tentArray[i] == tentArray[j]) {
                omp_unset_lock(&mylock1->lock);
                return;
            }
        } while (mylock1 != tentArray[j]);
        while (omp_test_lock(&(mylock2=tentArray[i])->lock)) {
            if (mylock2 == tentArray[i])
                goto locksOkay;
            omp_unset_lock(&mylock2->lock);
        }
        do {
            omp_unset_lock(&mylock1->lock);
            mylock1 = tentArray[i];
            omp_set_lock(&mylock1->lock);
            if (tentArray[i] == tentArray[j]) {
                omp_unset_lock(&mylock1->lock);
                return;
            }
        } while (mylock1 != tentArray[i]);
    }
    // now we have both locks, can't be identical!
    locksOkay:
    if (tentArray[i]->n < tentArray[j]->n) {
        integrate(tentArray[i], tentArray[j]);
        overwrite(tentArray[i], tentArray[j]);
    } else {
        integrate(tentArray[j], tentArray[i]);
        overwrite(tentArray[j], tentArray[i]);
    }
    //    printf("%2d UNLOCK %2d %2d\n", count, mylock1->minid, mylock2->minid);
    omp_unset_lock(&mylock1->lock);
    omp_unset_lock(&mylock2->lock);
}

void unifyTreeDbscan(int i, int j, int iCore, int jCore) {
    // merge if both are core objects or one is core object and the other is solitaire (n==1)
    if (tentArray[i] == tentArray[j])
        return;
    TreeEl * mylock1;
    TreeEl * mylock2;
    for(;;) {
        mylock1 = tentArray[i];
        omp_set_lock(&mylock1->lock);
        if (tentArray[i] == tentArray[j]) {
            omp_unset_lock(&mylock1->lock);
            return;
        }
        if (mylock1 == tentArray[i])
            break;
        omp_unset_lock(&mylock1->lock);
    }
    for (;;) {
        while (omp_test_lock(&(mylock2=tentArray[j])->lock)){
            if (mylock2 == tentArray[j])
                goto locksOkay;
            omp_unset_lock(&mylock2->lock);
        }
        do {
            omp_unset_lock(&mylock1->lock);
            mylock1 = tentArray[j];
            omp_set_lock(&mylock1->lock);
            if (tentArray[i] == tentArray[j]) {
                omp_unset_lock(&mylock1->lock);
                return;
            }
        } while (mylock1 != tentArray[j]);
        while (omp_test_lock(&(mylock2=tentArray[i])->lock)) {
            if (mylock2 == tentArray[i])
                goto locksOkay;
            omp_unset_lock(&mylock2->lock);
        }
        do {
            omp_unset_lock(&mylock1->lock);
            mylock1 = tentArray[i];
            omp_set_lock(&mylock1->lock);
            if (tentArray[i] == tentArray[j]) {
                omp_unset_lock(&mylock1->lock);
                return;
            }
        } while (mylock1 != tentArray[i]);
    }
    // now we have both locks, can't be identical!
    locksOkay:
    if (iCore && (jCore || tentArray[j]->n == 1) || jCore && tentArray[i]->n == 1) {
        if (tentArray[i]->n < tentArray[j]->n) {
            integrate(tentArray[i], tentArray[j]);
            overwrite(tentArray[i], tentArray[j]);
        } else {
            integrate(tentArray[j], tentArray[i]);
            overwrite(tentArray[j], tentArray[i]);
        }
    }
    //    printf("%2d UNLOCK %2d %2d\n", count, mylock1->minid, mylock2->minid);
    omp_unset_lock(&mylock1->lock);
    omp_unset_lock(&mylock2->lock);
}

void unifyTreeOrderedLock(int i, int j, int count) {
    //    printf("%2d LOCKni %2d %2d\n", count, i, tentArray[i]->minid);
    TreeEl * mylock1 = tentArray[i];
    TreeEl * mylock2 = tentArray[j];
    TreeEl * mylock3;
    for (;;) {
        if (mylock1 == mylock2)
            return;
        if (mylock1 < mylock2) {
            omp_set_lock(&mylock1->lock);
            if (tentArray[i] != mylock1 || (mylock2 = tentArray[j]) <= mylock1) {
                omp_unset_lock(&mylock1->lock);
                mylock1 = tentArray[i];
            } else {
                omp_set_lock(&mylock2->lock);
                while (mylock2 != tentArray[j] && (mylock3 = tentArray[j]) > mylock1) {
                    omp_unset_lock(&mylock2->lock);
                    mylock2 = mylock3;
                    omp_set_lock(&mylock2->lock);
                }
                if (mylock1 < mylock2)
                    break;
                omp_unset_lock(&mylock1->lock);
                omp_unset_lock(&mylock2->lock);
            }
        } else {
            omp_set_lock(&mylock2->lock);
            if (tentArray[j] != mylock2 || (mylock1 = tentArray[i]) <= mylock2) {
                omp_unset_lock(&mylock2->lock);
                mylock2 = tentArray[j];
            } else {
                omp_set_lock(&mylock1->lock);
                while (mylock1 != tentArray[i] && (mylock3 = tentArray[i]) > mylock2) {
                    omp_unset_lock(&mylock1->lock);
                    mylock1 = mylock3;
                    omp_set_lock(&mylock1->lock);
                }
                if (mylock1 > mylock2)
                    break;
                omp_unset_lock(&mylock1->lock);
                omp_unset_lock(&mylock2->lock);
            }
        }
    }
    if (tentArray[i]->n < tentArray[j]->n) {
        integrate(tentArray[i], tentArray[j]);
        overwrite(tentArray[i], tentArray[j]);
    } else {
        integrate(tentArray[j], tentArray[i]);
        overwrite(tentArray[j], tentArray[i]);
    }
    omp_unset_lock(&mylock1->lock);
    omp_unset_lock(&mylock2->lock);
}

void unifyTreeNoSync(int i, int j, int count) {
    if (tentArray[i] == tentArray[j])
        return;
    if (tentArray[i]->n < tentArray[j]->n) {
        integrate(tentArray[i], tentArray[j]);
        overwrite(tentArray[i], tentArray[j]);
    } else {
        integrate(tentArray[j], tentArray[i]);
        overwrite(tentArray[j], tentArray[i]);
    }
}
void unifyTreeSemiSync(int i, int j, int count) {
    if (tentArray[i] == tentArray[j])
        return;
    if (tentArray[i]->n < tentArray[j]->n) {
        omp_set_lock(&tentArray[j]->lock);
        integrate(tentArray[i], tentArray[j]);
        omp_unset_lock(&tentArray[j]->lock);
        overwrite(tentArray[i], tentArray[j]);
    } else {
        omp_set_lock(&tentArray[i]->lock);
        integrate(tentArray[j], tentArray[i]);
        omp_unset_lock(&tentArray[i]->lock);
        overwrite(tentArray[j], tentArray[i]);
    }
}
void printTree(int i){
    printf("%d: ",i);
    TreeEl *t= tentArray[i];
    while(t){
        printf("%d(%d) ", t->minid, t->n);
        t=t->left;
    }
    printf("\n");
}
void unifyTreeNoSync2(int i, int j, int count) {
//    int verbose = (count == 0 || count == 2 || count==95009 || count == 69327);//tentArray[i]->n == 2 && tentArray[j]->n == 2;
//    if (verbose) printf("*** ATTENTION ***\n%d %d %d\n", i,j,count);
//    if (verbose) printTree(i);
//    if (verbose) printTree(j);
    if (tentArray[i] == tentArray[j])
        return;
    if (tentArray[i]->n < tentArray[j]->n) {
        TreeEl *t = tentArray[i];
        TreeEl *s = t;
        while (t && t->left) {
            t = t->left;
//            if(verbose) printf("-> %d %ld ", t->minid, (long long) tentArray[j]);
            tentArray[t->minid] = tentArray[j];
        }
//        if(verbose)printf("\n");
//        if(verbose)printf("t: %d, tentArray[j]: %d\n", t->minid, tentArray[j]->minid);
        t->left = tentArray[j]->left;
        tentArray[j]->left = s;
        tentArray[j]->n += s->n;
        tentArray[s->minid] = tentArray[j];
        if (s->minid < tentArray[j]->minid) {
            int h = s->minid;
            s->minid = tentArray[j]->minid;
            tentArray[j]->minid = h;
        }
    } else {
        TreeEl *t = tentArray[j];
        TreeEl *s = t;
        while (t && t->left) {
            t = t->left;
//            if(verbose) printf("-> %d %ld ", t->minid, (long long) tentArray[j]);
            tentArray[t->minid] = tentArray[i];
        }
//        if(verbose)printf("\n");
//        if(verbose)printf("t: %d, tentArray[i]: %d\n", t->minid, tentArray[i]->minid);
        t->left = tentArray[i]->left;
//        if(verbose)printf("t->left: %d, tentArray[i]->left: %d\n", t->left->minid, tentArray[i]->left->minid);
        tentArray[i]->left = s;
        tentArray[i]->n += s->n;
        tentArray[s->minid] = tentArray[i];
        if (s->minid < tentArray[i]->minid) {
            int h = tentArray[i]->minid;
            tentArray[i]->minid = s->minid;
            s->minid = h;
        }
    }
//    if (verbose) printTree(i);
//    if (verbose) printTree(j);
//    if(tentArray[i]->minid > i || tentArray[j]->minid > j){
//        printf("Error %d %d %d\n", i,j,count);
//        exit(0);
//    }
//    if (verbose) exit(0);
}

void unifyList(int i, int j, int count) {
    if (listArray[i] == listArray[j])
        return;
    ListEl * mylock1;
    ListEl * mylock2;
    for(;;) {
        mylock1 = listArray[i];
        omp_set_lock(&mylock1->lock);
        if (listArray[i] == listArray[j]) {
            omp_unset_lock(&mylock1->lock);
            return;
        }
        if (mylock1 == listArray[i])
            break;
        omp_unset_lock(&mylock1->lock);
    }
    for (;;) {
        while (omp_test_lock(&(mylock2=listArray[j])->lock)){
            if (mylock2 == listArray[j])
                goto locksOkay;
            omp_unset_lock(&mylock2->lock);
        }
        do {
            omp_unset_lock(&mylock1->lock);
            mylock1 = listArray[j];
            omp_set_lock(&mylock1->lock);
            if (listArray[i] == listArray[j]) {
                omp_unset_lock(&mylock1->lock);
                return;
            }
        } while (mylock1 != listArray[j]);
        while (omp_test_lock(&(mylock2=listArray[i])->lock)) {
            if (mylock2 == listArray[i])
                goto locksOkay;
            omp_unset_lock(&mylock2->lock);
        }
        do {
            omp_unset_lock(&mylock1->lock);
            mylock1 = listArray[i];
            omp_set_lock(&mylock1->lock);
            if (listArray[i] == listArray[j]) {
                omp_unset_lock(&mylock1->lock);
                return;
            }
        } while (mylock1 != listArray[i]);
    }
    // now we have both locks, can't be identical!
    locksOkay:
    ListEl *t, *s, *r;
    if (listArray[i]->n < listArray[j]->n) {
        t = s = listArray[i];
        r = listArray[j];
    } else {
        t = s = listArray[j];
        r = listArray[i];
    }
    while (t && t->left) {
        t = t->left;
        listArray[t->minid] = r;
    }
    t->left = r->left;
    r->left = s;
    r->n += s->n;
    listArray[s->minid] = r;
    if (s->minid < r->minid) {
        int h = s->minid;
        s->minid = r->minid;
        r->minid = h;
    }
    omp_unset_lock(&mylock1->lock);
    omp_unset_lock(&mylock2->lock);
}

//long long unifyCounter;
//#define HISTO_SIZE 200
//int *histo;
void unifyListDbscan(int i, int j, int iCore, int jCore) {
    if (listArray[i] == listArray[j])
        return;
//    if(listArray[i] > listArray[j]) {
//        int h=i; i=j; j=h;
//        h=iCore; iCore=jCore; jCore=h;
//    }
    ListEl * mylock1;
    ListEl * mylock2;
    for(;;) {
        mylock1 = listArray[i];
        omp_set_lock(&mylock1->lock);
        if (listArray[i] == listArray[j]) {
            omp_unset_lock(&mylock1->lock);
            return;
        }
        if (mylock1 == listArray[i])
            break;
        omp_unset_lock(&mylock1->lock);
    }
    for (;;) {
        while (omp_test_lock(&(mylock2=listArray[j])->lock)){
            if (mylock2 == listArray[j])
                goto locksOkay;
            omp_unset_lock(&mylock2->lock);
        }
        do {
            omp_unset_lock(&mylock1->lock);
            mylock1 = listArray[j];
            omp_set_lock(&mylock1->lock);
            if (listArray[i] == listArray[j]) {
                omp_unset_lock(&mylock1->lock);
                return;
            }
        } while (mylock1 != listArray[j]);
        while (omp_test_lock(&(mylock2=listArray[i])->lock)) {
            if (mylock2 == listArray[i])
                goto locksOkay;
            omp_unset_lock(&mylock2->lock);
        }
        do {
            omp_unset_lock(&mylock1->lock);
            mylock1 = listArray[i];
            omp_set_lock(&mylock1->lock);
            if (listArray[i] == listArray[j]) {
                omp_unset_lock(&mylock1->lock);
                return;
            }
        } while (mylock1 != listArray[i]);
    }
    // now we have both locks, can't be identical!
    locksOkay:
    if (iCore && (jCore || listArray[j]->n == 1) || jCore && listArray[i]->n == 1) {
        ListEl *t, *s, *r;
        if (listArray[i]->n < listArray[j]->n) {
//            if(listArray[i]->n < HISTO_SIZE)
//                histo[listArray[i]->n]++;
//            else
//                histo[0]++;
            t = s = listArray[i];
            r = listArray[j];
        } else {
//            if(listArray[j]->n < HISTO_SIZE)
//                histo[listArray[j]->n]++;
//            else
//                histo[0]++;
            t = s = listArray[j];
            r = listArray[i];
        }
        while (t && t->left) {
//            unifyCounter++;
            t = t->left;
            listArray[t->minid] = r;
        }
        t->left = r->left;
        r->left = s;
        r->n += s->n;
        listArray[s->minid] = r;
        if (s->minid < r->minid) {
            int h = s->minid;
            s->minid = r->minid;
            r->minid = h;
        }
    }
    omp_unset_lock(&mylock1->lock);
    omp_unset_lock(&mylock2->lock);
}

void unifyListNoSync(int i, int j, int count) {
    if (listArray[i] == listArray[j])
        return;
    ListEl *t, *s, *r;
    if (listArray[i]->n < listArray[j]->n) {
        t = s = listArray[i];
        r = listArray[j];
    } else {
        t = s = listArray[j];
        r = listArray[i];
    }
    while (t && t->left) {
        t = t->left;
        listArray[t->minid] = r;
    }
    t->left = r->left;
    r->left = s;
    r->n += s->n;
    listArray[s->minid] = r;
    if (s->minid < r->minid) {
        int h = s->minid;
        s->minid = r->minid;
        r->minid = h;
    }
}

void unifyListCritical(int i, int j, int count) {
    if (listArray[i] == listArray[j])
        return;
#pragma omp critical
    {
        if (listArray[i] != listArray[j]) {
            ListEl *t, *s, *r;
            if (listArray[i]->n < listArray[j]->n) {
                t = s = listArray[i];
                r = listArray[j];
            } else {
                t = s = listArray[j];
                r = listArray[i];
            }
            while (t && t->left) {
                t = t->left;
                listArray[t->minid] = r;
            }
            t->left = r->left;
            r->left = s;
            r->n += s->n;
            listArray[s->minid] = r;
            if (s->minid < r->minid) {
                int h = s->minid;
                s->minid = r->minid;
                r->minid = h;
            }
        }
    }
}

void unifyTreeCritical(int i, int j, int count) {
    if (tentArray[i] == tentArray[j])
        return;
#pragma omp critical
    {
        if (tentArray[i] != tentArray[j])
            if (tentArray[i]->n < tentArray[j]->n) {
                integrate(tentArray[i], tentArray[j]);
                overwrite(tentArray[i], tentArray[j]);
            } else {
                integrate(tentArray[j], tentArray[i]);
                overwrite(tentArray[j], tentArray[i]);
            }
    }
}

void unifyTreeTest(int i, int j, int count) {
    if (tentArray[i] == tentArray[j])
        return;
    TreeEl * mylock1;
    TreeEl * mylock2;
    if (omp_test_lock(&(mylock1 = tentArray[i])->lock)) {
        if (omp_test_lock(&(mylock2 = tentArray[j])->lock)) {
            if (mylock1 == tentArray[i] && mylock2 == tentArray[j]) {
                if (tentArray[i]->n < tentArray[j]->n) {
                    integrate(tentArray[i], tentArray[j]);
                    overwrite(tentArray[i], tentArray[j]);
                } else {
                    integrate(tentArray[j], tentArray[i]);
                    overwrite(tentArray[j], tentArray[i]);
                }
            }
            omp_unset_lock(&mylock2->lock);
        }
        omp_unset_lock(&mylock1->lock);
    }
}

void testTree(int n, int m){
    tentArray=(TreeEl **) malloc(n * sizeof (TreeEl*));
    stop(); stopc();
    for (int i=0 ; i<n ; i++){
        tentArray[i]=(TreeEl *)malloc(sizeof(TreeEl));
        tentArray[i]->n=1;
        tentArray[i]->minid=i ;
        tentArray[i]->left = (TreeEl*)0;
        tentArray[i]->right = (TreeEl*)0;
        omp_init_lock(&tentArray[i]->lock);
    }
    int *testsched=(int*)malloc(2*m*sizeof(int));
    for (int i=0 ; i<m ; i++) {
        testsched[2*i]=drand48() * n;
        testsched[2*i+1]=drand48() * n;
    }
    printf("Setup     | %8.3f | %8.3f\n", stopc(), stop());
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++)
        for (int i = par * m / NUM_THREADS; i < (par + 1) * m / NUM_THREADS; i++) {
            int a1 = testsched[2 * i];
            int a2 = testsched[2 * i + 1];
            unifyTree(a1, a2, i);
        }
    for (int ii = 0; ii < n && ii < 100; ii += 10) {
        for (int i = ii; i < ii + 10 && i < n; i++)
            printf("%3d ", tentArray[i]->minid);
        for (int i = ii; i < ii + 10 && i < n; i++)
            if (omp_test_lock(&tentArray[i]->lock)) {
                printf("+ ");
                omp_unset_lock(&tentArray[i]->lock);
            } else printf("- ");
        printf("\n");
    }
    printf("\n");

    printf("Semaphor  | %8.3f | %8.3f\n", stopc(), stop());

    TreeEl ** tentArray2=(TreeEl **) malloc(n * sizeof (TreeEl*));
    memcpy(tentArray2, tentArray, n * sizeof (TreeEl*)) ;
    listArray = (ListEl **) malloc(n * sizeof (ListEl*));
    listArray[0] = (ListEl *) malloc(n * sizeof(ListEl));
    for (int i=0 ; i<n ; i++){
        listArray[i] = listArray[0] + i;
        listArray[i]->n=1;
        listArray[i]->minid=i ;
        listArray[i]->left = (ListEl*)0;
        omp_init_lock(&listArray[i]->lock);
    }
    printf("Setup     | %8.3f | %8.3f\n", stopc(), stop());
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++)
        for (int i = par * m / NUM_THREADS; i < (par + 1) * m / NUM_THREADS; i++) {
//    for(int i=0; i<m ; i++){
            int a1 = testsched[2 * i];
            int a2 = testsched[2 * i + 1];
            unifyList(a1, a2, i);
        }
    for (int ii=0; ii<n && ii<100 ; ii+=10){
    for (int i=ii ; i<ii+10 && i<n ; i++)
        printf("%3d ", listArray[i]->minid);
    for (int i=ii ; i<ii+10 && i<n ; i++)
        if(omp_test_lock(&listArray[i]->lock)){
            printf("+ ");
            omp_unset_lock(&listArray[i]->lock);
        } else printf("- ");
    printf("\n");
    }
    printf("\n");
    printf("Sequential| %8.3f | %8.3f\n", stopc(), stop());
    int count = 0;
    for (int i=0 ; i<n && count<20 ; i++)
        if(listArray[i]->minid != tentArray2[i]->minid)
            printf("%d. Error %d\n", ++count, i);
    printf("Number of Errors: %d\n", count);
}

int printout_vars(unsigned long long hilbert, int i, int j, int level, int clevel, int c, int lb0, int ub0, int lb1, int ub1) {
    if(clevel >= 10)
    printf("\t\t%ld %d %d %d %d %d %d %d %d %d\n", hilbert, i, j, level, clevel, c, lb0, ub0, lb1, ub1);
    return 1;
}

static inline void avxprint(vec a){
    double b[8];
    _mm512_storeu_pd(b,a);
    printf("                      %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
}
int nn;
int * L1;
int * L2;
int * L3;

int FGF_debug(int clevel, int lb0, int ub0, int lb1, int ub1){
    for(int i=lb0 ; i<ub0 ; i++)
        for (int j=lb1 ; j<ub1 ; j++)
            if(j>i && j < L1[i+nn/8]){
                if(!(ub1-1 > lb0 && lb1 < L1[((ub0-1)>>(clevel+1))+(nn>>(clevel+4))]))
                    printf("Shit... %d %d %d %d %d\n", clevel, lb0, ub0, lb1, ub1);
                return 1;
            }
    return 1;
}

#define transposeAVX512(r1,r2,r3,r4,r5,r6,r7,r8){\
            vec r1a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r1), mask1, _mm512_castpd_si512 (r2)));\
            vec r2a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r2), mask2, _mm512_castpd_si512 (r1)));\
            vec r3a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r3), mask1, _mm512_castpd_si512 (r4)));\
            vec r4a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r4), mask2, _mm512_castpd_si512 (r3)));\
            vec r5a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r5), mask1, _mm512_castpd_si512 (r6)));\
            vec r6a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r6), mask2, _mm512_castpd_si512 (r5)));\
            vec r7a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r7), mask1, _mm512_castpd_si512 (r8)));\
            vec r8a = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r8), mask2, _mm512_castpd_si512 (r7)));\
            vec r1b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r1a), mask3, _mm512_castpd_si512 (r3a)));\
            vec r3b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r3a), mask4, _mm512_castpd_si512 (r1a)));\
            vec r2b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r2a), mask3, _mm512_castpd_si512 (r4a)));\
            vec r4b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r4a), mask4, _mm512_castpd_si512 (r2a)));\
            vec r5b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r5a), mask3, _mm512_castpd_si512 (r7a)));\
            vec r7b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r7a), mask4, _mm512_castpd_si512 (r5a)));\
            vec r6b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r6a), mask3, _mm512_castpd_si512 (r8a)));\
            vec r8b = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r8a), mask4, _mm512_castpd_si512 (r6a)));\
            r1 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r1b), mask5, _mm512_castpd_si512 (r5b)));\
            r5 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r5b), mask6, _mm512_castpd_si512 (r1b)));\
            r2 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r2b), mask5, _mm512_castpd_si512 (r6b)));\
            r6 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r6b), mask6, _mm512_castpd_si512 (r2b)));\
            r3 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r3b), mask5, _mm512_castpd_si512 (r7b)));\
            r7 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r7b), mask6, _mm512_castpd_si512 (r3b)));\
            r4 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r4b), mask5, _mm512_castpd_si512 (r8b)));\
            r8 = _mm512_castsi512_pd (_mm512_permutex2var_epi64 (_mm512_castpd_si512 (r8b), mask6, _mm512_castpd_si512 (r4b)));\
            }
#define transposeAVX512i(r1,r2,r3,r4,r5,r6,r7,r8){\
            veci64 r1a = _mm512_permutex2var_epi64 (r1, mask1, r2);\
            veci64 r2a = _mm512_permutex2var_epi64 (r2, mask2, r1);\
            veci64 r3a = _mm512_permutex2var_epi64 (r3, mask1, r4);\
            veci64 r4a = _mm512_permutex2var_epi64 (r4, mask2, r3);\
            veci64 r5a = _mm512_permutex2var_epi64 (r5, mask1, r6);\
            veci64 r6a = _mm512_permutex2var_epi64 (r6, mask2, r5);\
            veci64 r7a = _mm512_permutex2var_epi64 (r7, mask1, r8);\
            veci64 r8a = _mm512_permutex2var_epi64 (r8, mask2, r7);\
            veci64 r1b = _mm512_permutex2var_epi64 (r1a, mask3, r3a);\
            veci64 r3b = _mm512_permutex2var_epi64 (r3a, mask4, r1a);\
            veci64 r2b = _mm512_permutex2var_epi64 (r2a, mask3, r4a);\
            veci64 r4b = _mm512_permutex2var_epi64 (r4a, mask4, r2a);\
            veci64 r5b = _mm512_permutex2var_epi64 (r5a, mask3, r7a);\
            veci64 r7b = _mm512_permutex2var_epi64 (r7a, mask4, r5a);\
            veci64 r6b = _mm512_permutex2var_epi64 (r6a, mask3, r8a);\
            veci64 r8b = _mm512_permutex2var_epi64 (r8a, mask4, r6a);\
            r1 = _mm512_permutex2var_epi64 (r1b, mask5, r5b);\
            r5 = _mm512_permutex2var_epi64 (r5b, mask6, r1b);\
            r2 = _mm512_permutex2var_epi64 (r2b, mask5, r6b);\
            r6 = _mm512_permutex2var_epi64 (r6b, mask6, r2b);\
            r3 = _mm512_permutex2var_epi64 (r3b, mask5, r7b);\
            r7 = _mm512_permutex2var_epi64 (r7b, mask6, r3b);\
            r4 = _mm512_permutex2var_epi64 (r4b, mask5, r8b);\
            r8 = _mm512_permutex2var_epi64 (r8b, mask6, r4b);\
            }

#define printAVX512(r){\
        _mm512_storeu_pd(partialresult, r);\
        printf("%6.2lf%6.2lf%6.2lf%6.2lf%6.2lf%6.2lf%6.2lf%6.2lf\n", \
                partialresult[0],partialresult[1],partialresult[2],partialresult[3],\
                partialresult[4],partialresult[5],partialresult[6],partialresult[7]);\
        }

int JoinCanonicAVX512(int n, int d, double epsilon, double* array) {
    long long result = 0;
    double epsilon2 = epsilon * epsilon;
    veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};
    veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};
    veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};
    veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};
    veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};
    veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};
    veci64 resultvec = _mm512_setzero_si512();
    //veci64 maskone = _mm512_set1_epi64(1ll);
    double partialresultX[72];
    double *partialresult = (double *) (((long long) partialresultX | 63ll) - 63ll);

    double *self = mallocA64(sizeof (double) * n);
    for (int i = 0; i < n / 8; i++) {
        register vec vi1 = _mm512_load_pd(array + d * 8 * i);
        register vec vi2 = _mm512_load_pd(array + d * (8 * i + 1));
        register vec vj = vi1;
        register vec sum1 = vi1 * vj;
        register vec sum9 = vi2 * vj;
        vj = vi2;
        register vec sum2 = vi1 * vj;
        register vec sum10 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 2));
        register vec sum3 = vi1 * vj;
        register vec sum11 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 3));
        register vec sum4 = vi1 * vj;
        register vec sum12 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 4));
        register vec sum5 = vi1 * vj;
        register vec sum13 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 5));
        register vec sum6 = vi1 * vj;
        register vec sum14 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 6));
        register vec sum7 = vi1 * vj;
        register vec sum15 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 7));
        register vec sum8 = vi1 * vj;
        register vec sum16 = vi2 * vj;
        for (register int k = 8; k < d; k += 8) {
            vi1 = _mm512_load_pd(array + d * 8 * i + k);
            vi2 = _mm512_load_pd(array + d * (8 * i + 1) + k);
            vj = vi1;
            sum1 += vi1 * vj;
            sum9 += vi2 * vj;
            vj = vi2;
            sum2 += vi1 * vj;
            sum10 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 2) + k);
            sum3 += vi1 * vj;
            sum11 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 3) + k);
            sum4 += vi1 * vj;
            sum12 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 4) + k);
            sum5 += vi1 * vj;
            sum13 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 5) + k);
            sum6 += vi1 * vj;
            sum14 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 6) + k);
            sum7 += vi1 * vj;
            sum15 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 7) + k);
            sum8 += vi1 * vj;
            sum16 += vi2 * vj;
        }
        //printAVX512(sum1);printAVX512(sum2);printAVX512(sum3);printAVX512(sum4);
        //printAVX512(sum5);printAVX512(sum6);printAVX512(sum7);printAVX512(sum8);
        transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
        //printf("\n");
        //printAVX512(sum1);printAVX512(sum2);printAVX512(sum3);printAVX512(sum4);
        //printAVX512(sum5);printAVX512(sum6);printAVX512(sum7);printAVX512(sum8);
        //exit(0);
        _mm512_store_pd(partialresult, sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8);
        transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
        _mm512_store_pd(partialresult + 8, sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16);
        vi1 = _mm512_load_pd(array + d * (8 * i + 2));
        vi2 = _mm512_load_pd(array + d * (8 * i + 3));
        vj = _mm512_load_pd(array + d * 8 * i);
        sum1 = vi1 * vj;
        sum9 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 1));
        sum2 = vi1 * vj;
        sum10 = vi2 * vj;
        vj = vi1;
        sum3 = vi1 * vj;
        sum11 = vi2 * vj;
        vj = vi2;
        sum4 = vi1 * vj;
        sum12 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 4));
        sum5 = vi1 * vj;
        sum13 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 5));
        sum6 = vi1 * vj;
        sum14 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 6));
        sum7 = vi1 * vj;
        sum15 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 7));
        sum8 = vi1 * vj;
        sum16 = vi2 * vj;
        for (register int k = 8; k < d; k += 8) {
            vi1 = _mm512_load_pd(array + d * (8 * i + 2) + k);
            vi2 = _mm512_load_pd(array + d * (8 * i + 3) + k);
            vj = _mm512_load_pd(array + d * 8 * i + k);
            sum1 += vi1 * vj;
            sum9 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 1) + k);
            sum2 += vi1 * vj;
            sum10 += vi2 * vj;
            vj = vi1;
            sum3 += vi1 * vj;
            sum11 += vi2 * vj;
            vj = vi2;
            sum4 += vi1 * vj;
            sum12 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 4) + k);
            sum5 += vi1 * vj;
            sum13 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 5) + k);
            sum6 += vi1 * vj;
            sum14 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 6) + k);
            sum7 += vi1 * vj;
            sum15 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 7) + k);
            sum8 += vi1 * vj;
            sum16 += vi2 * vj;
        }
        transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
        _mm512_store_pd(partialresult + 16, sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8);
        transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
        _mm512_store_pd(partialresult + 24, sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16);
        vi1 = _mm512_load_pd(array + d * (8 * i + 4));
        vi2 = _mm512_load_pd(array + d * (8 * i + 5));
        vj = _mm512_load_pd(array + d * 8 * i);
        sum1 = vi1 * vj;
        sum9 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 1));
        sum2 = vi1 * vj;
        sum10 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 2));
        sum3 = vi1 * vj;
        sum11 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 3));
        sum4 = vi1 * vj;
        sum12 = vi2 * vj;
        vj = vi1;
        sum5 = vi1 * vj;
        sum13 = vi2 * vj;
        vj = vi2;
        sum6 = vi1 * vj;
        sum14 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 6));
        sum7 = vi1 * vj;
        sum15 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 7));
        sum8 = vi1 * vj;
        sum16 = vi2 * vj;
        for (register int k = 8; k < d; k += 8) {
            vi1 = _mm512_load_pd(array + d * (8 * i + 4) + k);
            vi2 = _mm512_load_pd(array + d * (8 * i + 5) + k);
            vj = _mm512_load_pd(array + d * 8 * i + k);
            sum1 += vi1 * vj;
            sum9 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 1) + k);
            sum2 += vi1 * vj;
            sum10 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 2) + k);
            sum3 += vi1 * vj;
            sum11 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 3) + k);
            sum4 += vi1 * vj;
            sum12 += vi2 * vj;
            vj = vi1;
            sum5 += vi1 * vj;
            sum13 += vi2 * vj;
            vj = vi2;
            sum6 += vi1 * vj;
            sum14 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 6) + k);
            sum7 += vi1 * vj;
            sum15 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 7) + k);
            sum8 += vi1 * vj;
            sum16 += vi2 * vj;
        }
        transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
        _mm512_store_pd(partialresult + 32, sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8);
        transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
        _mm512_store_pd(partialresult + 40, sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16);
        vi1 = _mm512_load_pd(array + d * (8 * i + 6));
        vi2 = _mm512_load_pd(array + d * (8 * i + 7));
        vj = _mm512_load_pd(array + d * 8 * i);
        sum1 = vi1 * vj;
        sum9 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 1));
        sum2 = vi1 * vj;
        sum10 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 2));
        sum3 = vi1 * vj;
        sum11 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 3));
        sum4 = vi1 * vj;
        sum12 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 4));
        sum5 = vi1 * vj;
        sum13 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 5));
        sum6 = vi1 * vj;
        sum14 = vi2 * vj;
        vj = vi1;
        sum7 = vi1 * vj;
        sum15 = vi2 * vj;
        vj = vi2;
        sum8 = vi1 * vj;
        sum16 = vi2 * vj;
        for (register int k = 8; k < d; k += 8) {
            vi1 = _mm512_load_pd(array + d * (8 * i + 6) + k);
            vi2 = _mm512_load_pd(array + d * (8 * i + 7) + k);
            vj = _mm512_load_pd(array + d * 8 * i + k);
            sum1 += vi1 * vj;
            sum9 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 1) + k);
            sum2 += vi1 * vj;
            sum10 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 2) + k);
            sum3 += vi1 * vj;
            sum11 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 3) + k);
            sum4 += vi1 * vj;
            sum12 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 4) + k);
            sum5 += vi1 * vj;
            sum13 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 5) + k);
            sum6 += vi1 * vj;
            sum14 += vi2 * vj;
            vj = vi1;
            sum7 += vi1 * vj;
            sum15 += vi2 * vj;
            vj = vi2;
            sum8 += vi1 * vj;
            sum16 += vi2 * vj;
        }
        transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
        _mm512_store_pd(partialresult + 48, sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8);
        transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
        _mm512_store_pd(partialresult + 56, sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16);
        self[8 * i] = epsilon2 / 4 - partialresult[0] / 2;
        self[8 * i + 1] = epsilon2 / 4 - partialresult[9] / 2;
        self[8 * i + 2] = epsilon2 / 4 - partialresult[18] / 2;
        self[8 * i + 3] = epsilon2 / 4 - partialresult[27] / 2;
        self[8 * i + 4] = epsilon2 / 4 - partialresult[36] / 2;
        self[8 * i + 5] = epsilon2 / 4 - partialresult[45] / 2;
        self[8 * i + 6] = epsilon2 / 4 - partialresult[54] / 2;
        self[8 * i + 7] = epsilon2 / 4 - partialresult[63] / 2;
        if (partialresult[1] + self[8 * i] + self[8 * i + 1] >= 0) result++;
        if (partialresult[2] + self[8 * i] + self[8 * i + 2] >= 0) result++;
        if (partialresult[3] + self[8 * i] + self[8 * i + 3] >= 0) result++;
        if (partialresult[4] + self[8 * i] + self[8 * i + 4] >= 0) result++;
        if (partialresult[5] + self[8 * i] + self[8 * i + 5] >= 0) result++;
        if (partialresult[6] + self[8 * i] + self[8 * i + 6] >= 0) result++;
        if (partialresult[7] + self[8 * i] + self[8 * i + 7] >= 0) result++;
        if (partialresult[10] + self[8 * i + 1] + self[8 * i + 2] >= 0) result++;
        if (partialresult[11] + self[8 * i + 1] + self[8 * i + 3] >= 0) result++;
        if (partialresult[12] + self[8 * i + 1] + self[8 * i + 4] >= 0) result++;
        if (partialresult[13] + self[8 * i + 1] + self[8 * i + 5] >= 0) result++;
        if (partialresult[14] + self[8 * i + 1] + self[8 * i + 6] >= 0) result++;
        if (partialresult[15] + self[8 * i + 1] + self[8 * i + 7] >= 0) result++;
        if (partialresult[19] + self[8 * i + 2] + self[8 * i + 3] >= 0) result++;
        if (partialresult[20] + self[8 * i + 2] + self[8 * i + 4] >= 0) result++;
        if (partialresult[21] + self[8 * i + 2] + self[8 * i + 5] >= 0) result++;
        if (partialresult[22] + self[8 * i + 2] + self[8 * i + 6] >= 0) result++;
        if (partialresult[23] + self[8 * i + 2] + self[8 * i + 7] >= 0) result++;
        if (partialresult[28] + self[8 * i + 3] + self[8 * i + 4] >= 0) result++;
        if (partialresult[29] + self[8 * i + 3] + self[8 * i + 5] >= 0) result++;
        if (partialresult[30] + self[8 * i + 3] + self[8 * i + 6] >= 0) result++;
        if (partialresult[31] + self[8 * i + 3] + self[8 * i + 7] >= 0) result++;
        if (partialresult[37] + self[8 * i + 4] + self[8 * i + 5] >= 0) result++;
        if (partialresult[38] + self[8 * i + 4] + self[8 * i + 6] >= 0) result++;
        if (partialresult[39] + self[8 * i + 4] + self[8 * i + 7] >= 0) result++;
        if (partialresult[46] + self[8 * i + 5] + self[8 * i + 6] >= 0) result++;
        if (partialresult[47] + self[8 * i + 5] + self[8 * i + 7] >= 0) result++;
        if (partialresult[55] + self[8 * i + 6] + self[8 * i + 7] >= 0) result++;
    }
    printf("Result after self: %d\n", result);
    for (int i = 1; i < n / 8; i++)
        for (int j = 0; j < i; j++) {
            register vec vi1 = _mm512_load_pd(array + d * 8 * i);
            register vec vi2 = _mm512_load_pd(array + d * (8 * i + 1));
            register vec vj = _mm512_load_pd(array + d * 8 * j);
            register vec sum1 = vi1 * vj;
            register vec sum9 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 1));
            register vec sum2 = vi1 * vj;
            register vec sum10 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 2));
            register vec sum3 = vi1 * vj;
            register vec sum11 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 3));
            register vec sum4 = vi1 * vj;
            register vec sum12 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 4));
            register vec sum5 = vi1 * vj;
            register vec sum13 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 5));
            register vec sum6 = vi1 * vj;
            register vec sum14 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 6));
            register vec sum7 = vi1 * vj;
            register vec sum15 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 7));
            register vec sum8 = vi1 * vj;
            register vec sum16 = vi2 * vj;
            for (register int k = 8; k < d; k += 8) {
                vi1 = _mm512_load_pd(array + d * 8 * i + k);
                vi2 = _mm512_load_pd(array + d * (8 * i + 1) + k);
                vj = _mm512_load_pd(array + d * 8 * j + k);
                sum1 += vi1 * vj; // _mm512_fmadd_pd
                sum9 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                sum2 += vi1 * vj; // _mm512_fmadd_pd
                sum10 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                sum3 += vi1 * vj; // _mm512_fmadd_pd
                sum11 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                sum4 += vi1 * vj; // _mm512_fmadd_pd
                sum12 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                sum5 += vi1 * vj; // _mm512_fmadd_pd
                sum13 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                sum6 += vi1 * vj; // _mm512_fmadd_pd
                sum14 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                sum7 += vi1 * vj; // _mm512_fmadd_pd
                sum15 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                sum8 += vi1 * vj;
                sum16 += vi2 * vj;
            }
            transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
            vec hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                    _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i]);
            resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
            // _mm512_store_pd(partialresult, hhh) ;
            transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
            hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                    _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 1]);
            resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
            vi1 = _mm512_load_pd(array + d * (8 * i + 2));
            vi2 = _mm512_load_pd(array + d * (8 * i + 3));
            vj = _mm512_load_pd(array + d * 8 * j);
            sum1 = vi1 * vj;
            sum9 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 1));
            sum2 = vi1 * vj;
            sum10 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 2));
            sum3 = vi1 * vj;
            sum11 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 3));
            sum4 = vi1 * vj;
            sum12 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 4));
            sum5 = vi1 * vj;
            sum13 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 5));
            sum6 = vi1 * vj;
            sum14 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 6));
            sum7 = vi1 * vj;
            sum15 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 7));
            sum8 = vi1 * vj;
            sum16 = vi2 * vj;
            for (register int k = 8; k < d; k += 8) {
                vi1 = _mm512_load_pd(array + d * (8 * i + 2) + k);
                vi2 = _mm512_load_pd(array + d * (8 * i + 3) + k);
                vj = _mm512_load_pd(array + d * 8 * j + k);
                sum1 += vi1 * vj; // _mm512_fmadd_pd
                sum9 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                sum2 += vi1 * vj; // _mm512_fmadd_pd
                sum10 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                sum3 += vi1 * vj; // _mm512_fmadd_pd
                sum11 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                sum4 += vi1 * vj; // _mm512_fmadd_pd
                sum12 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                sum5 += vi1 * vj; // _mm512_fmadd_pd
                sum13 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                sum6 += vi1 * vj; // _mm512_fmadd_pd
                sum14 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                sum7 += vi1 * vj; // _mm512_fmadd_pd
                sum15 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                sum8 += vi1 * vj;
                sum16 += vi2 * vj;
            }
            transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
            hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                    _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 2]);
            resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
            // _mm512_store_pd(partialresult, hhh) ;
            transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
            hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                    _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 3]);
            resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
            vi1 = _mm512_load_pd(array + d * (8 * i + 4));
            vi2 = _mm512_load_pd(array + d * (8 * i + 5));
            vj = _mm512_load_pd(array + d * 8 * j);
            sum1 = vi1 * vj;
            sum9 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 1));
            sum2 = vi1 * vj;
            sum10 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 2));
            sum3 = vi1 * vj;
            sum11 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 3));
            sum4 = vi1 * vj;
            sum12 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 4));
            sum5 = vi1 * vj;
            sum13 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 5));
            sum6 = vi1 * vj;
            sum14 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 6));
            sum7 = vi1 * vj;
            sum15 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 7));
            sum8 = vi1 * vj;
            sum16 = vi2 * vj;
            for (register int k = 8; k < d; k += 8) {
                vi1 = _mm512_load_pd(array + d * (8 * i + 4) + k);
                vi2 = _mm512_load_pd(array + d * (8 * i + 5) + k);
                vj = _mm512_load_pd(array + d * 8 * j + k);
                sum1 += vi1 * vj; // _mm512_fmadd_pd
                sum9 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                sum2 += vi1 * vj; // _mm512_fmadd_pd
                sum10 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                sum3 += vi1 * vj; // _mm512_fmadd_pd
                sum11 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                sum4 += vi1 * vj; // _mm512_fmadd_pd
                sum12 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                sum5 += vi1 * vj; // _mm512_fmadd_pd
                sum13 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                sum6 += vi1 * vj; // _mm512_fmadd_pd
                sum14 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                sum7 += vi1 * vj; // _mm512_fmadd_pd
                sum15 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                sum8 += vi1 * vj;
                sum16 += vi2 * vj;
            }
            transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
            hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                    _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 4]);
            resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
            // _mm512_store_pd(partialresult, hhh) ;
            transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
            hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                    _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 5]);
            resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
            vi1 = _mm512_load_pd(array + d * (8 * i + 6));
            vi2 = _mm512_load_pd(array + d * (8 * i + 7));
            vj = _mm512_load_pd(array + d * 8 * j);
            sum1 = vi1 * vj;
            sum9 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 1));
            sum2 = vi1 * vj;
            sum10 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 2));
            sum3 = vi1 * vj;
            sum11 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 3));
            sum4 = vi1 * vj;
            sum12 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 4));
            sum5 = vi1 * vj;
            sum13 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 5));
            sum6 = vi1 * vj;
            sum14 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 6));
            sum7 = vi1 * vj;
            sum15 = vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * j + 7));
            sum8 = vi1 * vj;
            sum16 = vi2 * vj;
            for (register int k = 8; k < d; k += 8) {
                vi1 = _mm512_load_pd(array + d * (8 * i + 6) + k);
                vi2 = _mm512_load_pd(array + d * (8 * i + 7) + k);
                vj = _mm512_load_pd(array + d * 8 * j + k);
                sum1 += vi1 * vj; // _mm512_fmadd_pd
                sum9 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                sum2 += vi1 * vj; // _mm512_fmadd_pd
                sum10 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                sum3 += vi1 * vj; // _mm512_fmadd_pd
                sum11 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                sum4 += vi1 * vj; // _mm512_fmadd_pd
                sum12 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                sum5 += vi1 * vj; // _mm512_fmadd_pd
                sum13 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                sum6 += vi1 * vj; // _mm512_fmadd_pd
                sum14 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                sum7 += vi1 * vj; // _mm512_fmadd_pd
                sum15 += vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                sum8 += vi1 * vj;
                sum16 += vi2 * vj;
            }
            transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
            hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                    _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 6]);
            resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
            // _mm512_store_pd(partialresult, hhh) ;
            transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
            hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                    _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 7]);
            resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
            result += 64;
        }

    long long resultvecX[8];
    _mm512_storeu_si512((void *) resultvecX, resultvec);
    //    printf("%ld %ld %ld %ld %ld %ld %ld %ld\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
    //    printf("%lx %lx %lx %lx %lx %lx %lx %lx\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
    result -= resultvecX[0] + resultvecX[1] + resultvecX[2] + resultvecX[3] + resultvecX[4] + resultvecX[5] + resultvecX[6] + resultvecX[7];
    printf("Result after all: %ld\n", 2 * result + n);
    return 2 * result + n;
}
double EGO_epsilon;
int EGO_d;
int min_size_qsort;

#define swap(i,j) {memcpy(b, (i), size); memcpy((i), (j), size); memcpy((j), b, size);}

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

int num_def_entries;
void ** start_def;
int * num_def;

void def_qsort_rek (void* l, size_t num, size_t size, int (*compar)(const void*,const void*)){
    printf ("Qsort %ld %d %d %d\n", (long long)l, num, omp_get_num_threads(), omp_get_thread_num());
    if(num <= min_size_qsort){
        if(num_def_entries > 0 && num < min_size_qsort/10){
            num_def[num_def_entries-1] += num+1 ;
        } else {
            start_def[num_def_entries] = l;
            num_def[num_def_entries++] = num;
        }
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
    def_qsort_rek (l, ((long long)i-(long long)l) / size, size, compar);
    def_qsort_rek (i+size, num - ((long long)i-(long long)l)/size - 1, size, compar);
}

void def_qsort (void* l, size_t num, size_t size, int (*compar)(const void*,const void*)){
    num_def_entries = 0;
    start_def = (void **) malloc (NUM_THREADS * 10 * sizeof(void *));
    num_def = (int *) malloc (NUM_THREADS * 10 * sizeof(int));
    def_qsort_rek (l, num, size, compar) ;
    printf("tasks: %d (%d)\n", num_def_entries, num_def[0]) ;
    printf("deferred q | %8.3f | %8.3f\n", stopc(), stop());

#pragma omp parallel for
    for(int i=0 ; i<num_def_entries ; i++)
        qsort(start_def[i], num_def[i], size, compar) ;
}

int epsilonGridCompare(const void *a, const void *b) {
    double * A = (double *) a;
    double * B = (double *) b;
    for (int i = 0; i < EGO_d; i++) {
        int h = (int) (floor(A[i] / EGO_epsilon) - floor(B[i] / EGO_epsilon));
        if (h != 0) return h;
    }
    return 0;
}

void epsilonGridOrdering(int n, int d, double epsilon, double* array) {
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

void epsilonGridFillList3(int n, int d, double epsilon, double * array, int *L1, int *L2, int *L3) {
    int a1, a2, a3;
    a1 = a2 = a3 = 0;
    int back1 = 0, back2 = 0, back3 = 0;
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
    printf("Backsteps: %d %d %d\n", back1, back2, back3);
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
#define max(a,b) ((a)>(b) ? (a) : (b))
#define min(a,b) ((a)<(b) ? (a) : (b))

void epsilonGridCompleteListMax(int n, int *list) {
    int m = ceilpowtwo(n);
    for (int i = n; i < m; i++)
        list[i+m] = (n+7)/8;
    if (n % 2)
        list[n] = list[2 * n - 1];
    for (int i = n - 1; i > 0; i--)
        list[i] = max(list[2 * i], list[2 * i + 1]);
}

void epsilonGridCompleteListMin(int n, int *list) {
    int m = ceilpowtwo(n);
    for (int i = n; i < m; i++)
        list[i+m] = 0;
    if (n % 2)
        list[n] = list[2 * n - 1];
    for (int i = n - 1; i > 0; i--)
        list[i] = min(list[2 * i], list[2 * i + 1]);
}

int joinCanonicAVX1(int n, int d, double epsilon, double* array) {
    int result = 0;
    double partialresultA[24];
    double *partialresult = (double*) (((unsigned long long) partialresultA | 63ull) + 1);
    double epsilon2 = epsilon*epsilon;
    stopc();
    stop();
    double *self = (double *) malloc(sizeof (double) * n);
    for (int i = 0; i < n / 4; i++) {
        register vec4 vi1 = _mm256_load_pd(array + d * 4 * i);
        register vec4 vi2 = _mm256_load_pd(array + d * (4 * i + 1));
        register vec4 vj = vi1;
        register vec4 sum1 = vi1 * vj;
        register vec4 sum5 = vi2 * vj;
        vj = vi2;
        register vec4 sum2 = vi1 * vj;
        register vec4 sum6 = vi2 * vj;
        vj = _mm256_load_pd(array + d * (4 * i + 2));
        register vec4 sum3 = vi1 * vj;
        register vec4 sum7 = vi2 * vj;
        vj = _mm256_load_pd(array + d * (4 * i + 3));
        register vec4 sum4 = vi1 * vj;
        register vec4 sum8 = vi2 * vj;
        for (register int k = 4; k < d; k += 4) {
            vi1 = _mm256_load_pd(array + d * 4 * i + k);
            vi2 = _mm256_load_pd(array + d * (4 * i + 1) + k);
            vj = vi1;
            sum1 += vi1 * vj;
            sum5 += vi2 * vj;
            vj = vi2;
            sum2 += vi1 * vj;
            sum6 += vi2 * vj;
            vj = _mm256_load_pd(array + d * (4 * i + 2) + k);
            sum3 += vi1 * vj;
            sum7 += vi2 * vj;
            vj = _mm256_load_pd(array + d * (4 * i + 3) + k);
            sum4 += vi1 * vj;
            sum8 += vi2 * vj;
        }
        vec4 sumab = _mm256_hadd_pd(sum1, sum2);
        vec4 sumcd = _mm256_hadd_pd(sum3, sum4);
        sumab = _mm256_blend_pd(sumab, sumcd, 0b1100)
                + _mm256_permute2f128_pd(sumab, sumcd, 0x21);
        //sumab += _mm256_load_pd(self + j * 4);
        //sumab += _mm256_broadcast_sd(self + i * 4);
        // now the sign bit is 0 iff (i*4, j*4+{0,1,2,3}) is join result
        _mm256_store_pd(partialresult, sumab);
        sumab = _mm256_hadd_pd(sum5, sum6);
        sumcd = _mm256_hadd_pd(sum7, sum8);
        sumab = _mm256_blend_pd(sumab, sumcd, 0b1100)
                + _mm256_permute2f128_pd(sumab, sumcd, 0x21);
        //sumab += _mm256_load_pd(self + j * 4);
        //sumab += _mm256_broadcast_sd(self + (i * 4 + 1));
        // now the sign bit is 0 iff (i*4+1, j*4+{0,1,2,3}) is join result
        _mm256_store_pd(partialresult + 4, sumab);
        vi1 = _mm256_load_pd(array + d * (4 * i + 2));
        vi2 = _mm256_load_pd(array + d * (4 * i + 3));
        vj = _mm256_load_pd(array + d * 4 * i);
        sum1 = vi1 * vj;
        sum5 = vi2 * vj;
        vj = _mm256_load_pd(array + d * (4 * i + 1));
        sum2 = vi1 * vj;
        sum6 = vi2 * vj;
        vj = vi1;
        sum3 = vi1 * vj;
        sum7 = vi2 * vj;
        vj = vi2;
        sum4 = vi1 * vj;
        sum8 = vi2 * vj;
        for (register int k = 4; k < d; k += 4) {
            vi1 = _mm256_load_pd(array + d * (4 * i + 2) + k);
            vi2 = _mm256_load_pd(array + d * (4 * i + 3) + k);
            vj = _mm256_load_pd(array + d * 4 * i + k);
            sum1 += vi1 * vj;
            sum5 += vi2 * vj;
            vj = _mm256_load_pd(array + d * (4 * i + 1) + k);
            sum2 += vi1 * vj;
            sum6 += vi2 * vj;
            vj = vi1;
            sum3 += vi1 * vj;
            sum7 += vi2 * vj;
            vj = vi2;
            sum4 += vi1 * vj;
            sum8 += vi2 * vj;
        }
        sumab = _mm256_hadd_pd(sum1, sum2);
        sumcd = _mm256_hadd_pd(sum3, sum4);
        sumab = _mm256_blend_pd(sumab, sumcd, 0b1100)
                + _mm256_permute2f128_pd(sumab, sumcd, 0x21);
        //sumab += _mm256_load_pd(self + j * 4);
        //sumab += _mm256_broadcast_sd(self + (4 * i + 2));
        // now the sign bit is 0 iff (i*4+2, j*4+{0,1,2,3}) is join result
        _mm256_store_pd(partialresult + 8, sumab);
        sumab = _mm256_hadd_pd(sum5, sum6);
        sumcd = _mm256_hadd_pd(sum7, sum8);
        sumab = _mm256_blend_pd(sumab, sumcd, 0b1100)
                + _mm256_permute2f128_pd(sumab, sumcd, 0x21);
        //sumab += _mm256_load_pd(self + j * 4);
        //sumab += _mm256_broadcast_sd(self + (4 * i + 3));
        // now the sign bit is 0 iff (i*4+3, j*4+{0,1,2,3}) is join result
        _mm256_store_pd(partialresult + 12, sumab);
        // the following could be done with _mm512_ror_epi(sumab, 63); _mm512_xor(sumab, 1); _mm512_add_epi(sumab, xxx)
        self[4 * i] = epsilon2 / 4 - partialresult[0] / 2;
        self[4 * i + 1] = epsilon2 / 4 - partialresult[5] / 2;
        self[4 * i + 2] = epsilon2 / 4 - partialresult[10] / 2;
        self[4 * i + 3] = epsilon2 / 4 - partialresult[15] / 2;
        if (partialresult[1] + self[4 * i] + self[4 * i + 1] >= 0) result++;
        if (partialresult[2] + self[4 * i] + self[4 * i + 2] >= 0) result++;
        if (partialresult[3] + self[4 * i] + self[4 * i + 3] >= 0) result++;
        if (partialresult[6] + self[4 * i + 1] + self[4 * i + 2] >= 0) result++;
        if (partialresult[7] + self[4 * i + 1] + self[4 * i + 3] >= 0) result++;
        if (partialresult[11] + self[4 * i + 2] + self[4 * i + 3] >= 0) result++;
    }
    printf("Self       | %8.3f | %8.3f\n", stopc(), stop());

    // Determine Euclidean distances between diagonal (excluded) and L1
    int squarecand = 0;
    int filteredout = 0;
    unsigned long long lasthilbert = 2ull;
    int samecounter = 0;
    //    int i=0; int j=0; FGF_HILBERT_FOR(i, j, n/4, n/4, j>i && j < L1[i+nn/4], /*printout_vars(HILLOOP_hilbert, i, j, HILLOOP_level, FURHIL_clevel, HILLOOP_c, FURHIL_lb0, FURHIL_ub0, FURHIL_lb1, FURHIL_ub1) &&*/ FURHIL_ub1-1 > FURHIL_lb0 && FURHIL_lb1 < L1[((FURHIL_ub0-1)>>(FURHIL_clevel+1))+(nn>>(FURHIL_clevel+3))]){
    for (int i = 0; i < n / 4; i++) for (int j = i + 1; j < n / 4; j++) {
            squarecand++;
            register vec4 vi1 = _mm256_load_pd(array + d * 4 * i);
            register vec4 vi2 = _mm256_load_pd(array + d * (4 * i + 1));
            register vec4 vj = _mm256_load_pd(array + d * 4 * j);
            register vec4 sum1 = vi1 * vj;
            register vec4 sum5 = vi2 * vj;
            vj = _mm256_load_pd(array + d * (4 * j + 1));
            register vec4 sum2 = vi1 * vj;
            register vec4 sum6 = vi2 * vj;
            vj = _mm256_load_pd(array + d * (4 * j + 2));
            register vec4 sum3 = vi1 * vj;
            register vec4 sum7 = vi2 * vj;
            vj = _mm256_load_pd(array + d * (4 * j + 3));
            register vec4 sum4 = vi1 * vj;
            register vec4 sum8 = vi2 * vj;
            for (register int k = 4; k < d; k += 4) {
                vi1 = _mm256_load_pd(array + d * 4 * i + k);
                vi2 = _mm256_load_pd(array + d * (4 * i + 1) + k);
                vj = _mm256_load_pd(array + d * 4 * j + k);
                sum1 += vi1 * vj;
                sum5 += vi2 * vj;
                vj = _mm256_load_pd(array + d * (4 * j + 1) + k);
                sum2 += vi1 * vj;
                sum6 += vi2 * vj;
                vj = _mm256_load_pd(array + d * (4 * j + 2) + k);
                sum3 += vi1 * vj;
                sum7 += vi2 * vj;
                vj = _mm256_load_pd(array + d * (4 * j + 3) + k);
                sum4 += vi1 * vj;
                sum8 += vi2 * vj;
            }
            vec4 sumab = _mm256_hadd_pd(sum1, sum2);
            vec4 sumcd = _mm256_hadd_pd(sum3, sum4);
            sumab = _mm256_blend_pd(sumab, sumcd, 0b1100)
                    + _mm256_permute2f128_pd(sumab, sumcd, 0x21);
            sumab += _mm256_load_pd(self + j * 4);
            sumab += _mm256_broadcast_sd(self + i * 4);
            // now the sign bit is 0 iff (i*4, j*4+{0,1,2,3}) is join result
            _mm256_store_pd(partialresult, sumab);
            sumab = _mm256_hadd_pd(sum5, sum6);
            sumcd = _mm256_hadd_pd(sum7, sum8);
            sumab = _mm256_blend_pd(sumab, sumcd, 0b1100)
                    + _mm256_permute2f128_pd(sumab, sumcd, 0x21);
            sumab += _mm256_load_pd(self + j * 4);
            sumab += _mm256_broadcast_sd(self + (i * 4 + 1));
            // now the sign bit is 0 iff (i*4+1, j*4+{0,1,2,3}) is join result
            _mm256_store_pd(partialresult + 4, sumab);
            vi1 = _mm256_load_pd(array + d * (4 * i + 2));
            vi2 = _mm256_load_pd(array + d * (4 * i + 3));
            vj = _mm256_load_pd(array + d * 4 * j);
            sum1 = vi1 * vj;
            sum5 = vi2 * vj;
            vj = _mm256_load_pd(array + d * (4 * j + 1));
            sum2 = vi1 * vj;
            sum6 = vi2 * vj;
            vj = _mm256_load_pd(array + d * (4 * j + 2));
            sum3 = vi1 * vj;
            sum7 = vi2 * vj;
            vj = _mm256_load_pd(array + d * (4 * j + 3));
            sum4 = vi1 * vj;
            sum8 = vi2 * vj;
            for (register int k = 4; k < d; k += 4) {
                vi1 = _mm256_load_pd(array + d * (4 * i + 2) + k);
                vi2 = _mm256_load_pd(array + d * (4 * i + 3) + k);
                vj = _mm256_load_pd(array + d * 4 * j + k);
                sum1 += vi1 * vj;
                sum5 += vi2 * vj;
                vj = _mm256_load_pd(array + d * (4 * j + 1) + k);
                sum2 += vi1 * vj;
                sum6 += vi2 * vj;
                vj = _mm256_load_pd(array + d * (4 * j + 2) + k);
                sum3 += vi1 * vj;
                sum7 += vi2 * vj;
                vj = _mm256_load_pd(array + d * (4 * j + 3) + k);
                sum4 += vi1 * vj;
                sum8 += vi2 * vj;
            }
            sumab = _mm256_hadd_pd(sum1, sum2);
            sumcd = _mm256_hadd_pd(sum3, sum4);
            sumab = _mm256_blend_pd(sumab, sumcd, 0b1100)
                    + _mm256_permute2f128_pd(sumab, sumcd, 0x21);
            sumab += _mm256_load_pd(self + j * 4);
            sumab += _mm256_broadcast_sd(self + (4 * i + 2));
            // now the sign bit is 0 iff (i*4+2, j*4+{0,1,2,3}) is join result
            _mm256_store_pd(partialresult + 8, sumab);
            sumab = _mm256_hadd_pd(sum5, sum6);
            sumcd = _mm256_hadd_pd(sum7, sum8);
            sumab = _mm256_blend_pd(sumab, sumcd, 0b1100)
                    + _mm256_permute2f128_pd(sumab, sumcd, 0x21);
            sumab += _mm256_load_pd(self + j * 4);
            sumab += _mm256_broadcast_sd(self + (4 * i + 3));
            // now the sign bit is 0 iff (i*4+3, j*4+{0,1,2,3}) is join result
            _mm256_store_pd(partialresult + 12, sumab);
            // the following could be done with _mm512_ror_epi(sumab, 63); _mm512_xor(sumab, 1); _mm512_add_epi(sumab, xxx)
            for (int k = 0; k < 16; k++)
                if (partialresult[k] >= 0) {
                    result++;
                }
            //} else {filteredout++; if (lasthilbert == HILLOOP_hilbert) {samecounter ++; if (samecounter>=64||HILLOOP_hilbert == 50183)printf ("error %d %ld %d %d L1: %d %d %d %d %d %d\n", samecounter, HILLOOP_hilbert, i, j, L1[i+nn/4], L1[i/2+nn/8], L1[i/4+nn/16], L1[i/8+nn/32], L1[i/16+nn/64], L1[i/32+nn/128]);}else{/*printf("%d %ld\n", samecounter, lasthilbert); */lasthilbert = HILLOOP_hilbert; samecounter =1;}
        } // FGF_HILBERT_END(i,j);
    printf("4x4 Square Candidates: %d (%d)\n\n", squarecand, filteredout);
    printf("Join AVX   | %8.3f | %8.3f\n", stopc(), stop());
    return 2 * result + n;
}

int epsilonGridJoinCanonicAVX512(int n, int d, double epsilon, double* array) {
    long long result = 0;
    double epsilon2 = epsilon * epsilon;
    veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};
    veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};
    veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};
    veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};
    veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};
    veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};
    //veci64 maskone = _mm512_set1_epi64(1ll);
    double partialresultX[72];
    double *partialresult = (double *) (((long long) partialresultX | 63ll) - 63ll);

    double *self = mallocA64(sizeof (double) * n);
    stopc();
    stop();
    epsilonGridOrdering(n, d, epsilon, array);
    printf("Sort       | %8.3f | %8.3f\n", stopc(), stop());
    /*int*/ nn = ceilpowtwo(n);
    printf("nn=%d\n",nn);
    /*int * */ L1 = (int *) mallocA64(sizeof (int) * 2 * nn);
    /*int * */ L2 = (int *) mallocA64(sizeof (int) * 2 * nn);
    /*int * */ L3 = (int *) mallocA64(sizeof (int) * 2 * nn);
    epsilonGridFillList3(n, d, epsilon, array, L1 + nn, L2 + nn, L3 + nn);
    for (int i = 0; i < n; i++)
        L1[nn + i] = (L1[nn + i] + 7) / 8;
    for (int i = 0; i < n; i++)
        L2[nn + i] = L2[nn + i] / 8;
    for (int i = 0; i < n; i++)
        L3[nn + i] = (L3[nn + i] + 7) / 8;
    epsilonGridCompleteListMax(nn, L1);
    epsilonGridCompleteListMin(nn, L2);
    epsilonGridCompleteListMax(nn, L3);
    printf("FillList   | %8.3f | %8.3f\n", stopc(), stop());
    for (int i = 0; i < n / 8; i++) {
        register vec vi1 = _mm512_load_pd(array + d * 8 * i);
        register vec vi2 = _mm512_load_pd(array + d * (8 * i + 1));
        register vec vj = vi1;
        register vec sum1 = vi1 * vj;
        register vec sum9 = vi2 * vj;
        vj = vi2;
        register vec sum2 = vi1 * vj;
        register vec sum10 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 2));
        register vec sum3 = vi1 * vj;
        register vec sum11 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 3));
        register vec sum4 = vi1 * vj;
        register vec sum12 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 4));
        register vec sum5 = vi1 * vj;
        register vec sum13 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 5));
        register vec sum6 = vi1 * vj;
        register vec sum14 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 6));
        register vec sum7 = vi1 * vj;
        register vec sum15 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 7));
        register vec sum8 = vi1 * vj;
        register vec sum16 = vi2 * vj;
        for (register int k = 8; k < d; k += 8) {
            vi1 = _mm512_load_pd(array + d * 8 * i + k);
            vi2 = _mm512_load_pd(array + d * (8 * i + 1) + k);
            vj = vi1;
            sum1 += vi1 * vj;
            sum9 += vi2 * vj;
            vj = vi2;
            sum2 += vi1 * vj;
            sum10 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 2) + k);
            sum3 += vi1 * vj;
            sum11 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 3) + k);
            sum4 += vi1 * vj;
            sum12 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 4) + k);
            sum5 += vi1 * vj;
            sum13 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 5) + k);
            sum6 += vi1 * vj;
            sum14 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 6) + k);
            sum7 += vi1 * vj;
            sum15 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 7) + k);
            sum8 += vi1 * vj;
            sum16 += vi2 * vj;
        }
        //printAVX512(sum1);printAVX512(sum2);printAVX512(sum3);printAVX512(sum4);
        //printAVX512(sum5);printAVX512(sum6);printAVX512(sum7);printAVX512(sum8);
        transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
        //printf("\n");
        //printAVX512(sum1);printAVX512(sum2);printAVX512(sum3);printAVX512(sum4);
        //printAVX512(sum5);printAVX512(sum6);printAVX512(sum7);printAVX512(sum8);
        //exit(0);
        _mm512_store_pd(partialresult, sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8);
        transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
        _mm512_store_pd(partialresult + 8, sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16);
        vi1 = _mm512_load_pd(array + d * (8 * i + 2));
        vi2 = _mm512_load_pd(array + d * (8 * i + 3));
        vj = _mm512_load_pd(array + d * 8 * i);
        sum1 = vi1 * vj;
        sum9 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 1));
        sum2 = vi1 * vj;
        sum10 = vi2 * vj;
        vj = vi1;
        sum3 = vi1 * vj;
        sum11 = vi2 * vj;
        vj = vi2;
        sum4 = vi1 * vj;
        sum12 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 4));
        sum5 = vi1 * vj;
        sum13 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 5));
        sum6 = vi1 * vj;
        sum14 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 6));
        sum7 = vi1 * vj;
        sum15 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 7));
        sum8 = vi1 * vj;
        sum16 = vi2 * vj;
        for (register int k = 8; k < d; k += 8) {
            vi1 = _mm512_load_pd(array + d * (8 * i + 2) + k);
            vi2 = _mm512_load_pd(array + d * (8 * i + 3) + k);
            vj = _mm512_load_pd(array + d * 8 * i + k);
            sum1 += vi1 * vj;
            sum9 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 1) + k);
            sum2 += vi1 * vj;
            sum10 += vi2 * vj;
            vj = vi1;
            sum3 += vi1 * vj;
            sum11 += vi2 * vj;
            vj = vi2;
            sum4 += vi1 * vj;
            sum12 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 4) + k);
            sum5 += vi1 * vj;
            sum13 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 5) + k);
            sum6 += vi1 * vj;
            sum14 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 6) + k);
            sum7 += vi1 * vj;
            sum15 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 7) + k);
            sum8 += vi1 * vj;
            sum16 += vi2 * vj;
        }
        transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
        _mm512_store_pd(partialresult + 16, sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8);
        transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
        _mm512_store_pd(partialresult + 24, sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16);
        vi1 = _mm512_load_pd(array + d * (8 * i + 4));
        vi2 = _mm512_load_pd(array + d * (8 * i + 5));
        vj = _mm512_load_pd(array + d * 8 * i);
        sum1 = vi1 * vj;
        sum9 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 1));
        sum2 = vi1 * vj;
        sum10 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 2));
        sum3 = vi1 * vj;
        sum11 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 3));
        sum4 = vi1 * vj;
        sum12 = vi2 * vj;
        vj = vi1;
        sum5 = vi1 * vj;
        sum13 = vi2 * vj;
        vj = vi2;
        sum6 = vi1 * vj;
        sum14 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 6));
        sum7 = vi1 * vj;
        sum15 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 7));
        sum8 = vi1 * vj;
        sum16 = vi2 * vj;
        for (register int k = 8; k < d; k += 8) {
            vi1 = _mm512_load_pd(array + d * (8 * i + 4) + k);
            vi2 = _mm512_load_pd(array + d * (8 * i + 5) + k);
            vj = _mm512_load_pd(array + d * 8 * i + k);
            sum1 += vi1 * vj;
            sum9 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 1) + k);
            sum2 += vi1 * vj;
            sum10 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 2) + k);
            sum3 += vi1 * vj;
            sum11 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 3) + k);
            sum4 += vi1 * vj;
            sum12 += vi2 * vj;
            vj = vi1;
            sum5 += vi1 * vj;
            sum13 += vi2 * vj;
            vj = vi2;
            sum6 += vi1 * vj;
            sum14 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 6) + k);
            sum7 += vi1 * vj;
            sum15 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 7) + k);
            sum8 += vi1 * vj;
            sum16 += vi2 * vj;
        }
        transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
        _mm512_store_pd(partialresult + 32, sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8);
        transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
        _mm512_store_pd(partialresult + 40, sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16);
        vi1 = _mm512_load_pd(array + d * (8 * i + 6));
        vi2 = _mm512_load_pd(array + d * (8 * i + 7));
        vj = _mm512_load_pd(array + d * 8 * i);
        sum1 = vi1 * vj;
        sum9 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 1));
        sum2 = vi1 * vj;
        sum10 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 2));
        sum3 = vi1 * vj;
        sum11 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 3));
        sum4 = vi1 * vj;
        sum12 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 4));
        sum5 = vi1 * vj;
        sum13 = vi2 * vj;
        vj = _mm512_load_pd(array + d * (8 * i + 5));
        sum6 = vi1 * vj;
        sum14 = vi2 * vj;
        vj = vi1;
        sum7 = vi1 * vj;
        sum15 = vi2 * vj;
        vj = vi2;
        sum8 = vi1 * vj;
        sum16 = vi2 * vj;
        for (register int k = 8; k < d; k += 8) {
            vi1 = _mm512_load_pd(array + d * (8 * i + 6) + k);
            vi2 = _mm512_load_pd(array + d * (8 * i + 7) + k);
            vj = _mm512_load_pd(array + d * 8 * i + k);
            sum1 += vi1 * vj;
            sum9 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 1) + k);
            sum2 += vi1 * vj;
            sum10 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 2) + k);
            sum3 += vi1 * vj;
            sum11 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 3) + k);
            sum4 += vi1 * vj;
            sum12 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 4) + k);
            sum5 += vi1 * vj;
            sum13 += vi2 * vj;
            vj = _mm512_load_pd(array + d * (8 * i + 5) + k);
            sum6 += vi1 * vj;
            sum14 += vi2 * vj;
            vj = vi1;
            sum7 += vi1 * vj;
            sum15 += vi2 * vj;
            vj = vi2;
            sum8 += vi1 * vj;
            sum16 += vi2 * vj;
        }
        transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
        _mm512_store_pd(partialresult + 48, sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8);
        transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
        _mm512_store_pd(partialresult + 56, sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16);
        self[8 * i] = epsilon2 / 4 - partialresult[0] / 2;
        self[8 * i + 1] = epsilon2 / 4 - partialresult[9] / 2;
        self[8 * i + 2] = epsilon2 / 4 - partialresult[18] / 2;
        self[8 * i + 3] = epsilon2 / 4 - partialresult[27] / 2;
        self[8 * i + 4] = epsilon2 / 4 - partialresult[36] / 2;
        self[8 * i + 5] = epsilon2 / 4 - partialresult[45] / 2;
        self[8 * i + 6] = epsilon2 / 4 - partialresult[54] / 2;
        self[8 * i + 7] = epsilon2 / 4 - partialresult[63] / 2;
        if (partialresult[1] + self[8 * i] + self[8 * i + 1] >= 0) result++;
        if (partialresult[2] + self[8 * i] + self[8 * i + 2] >= 0) result++;
        if (partialresult[3] + self[8 * i] + self[8 * i + 3] >= 0) result++;
        if (partialresult[4] + self[8 * i] + self[8 * i + 4] >= 0) result++;
        if (partialresult[5] + self[8 * i] + self[8 * i + 5] >= 0) result++;
        if (partialresult[6] + self[8 * i] + self[8 * i + 6] >= 0) result++;
        if (partialresult[7] + self[8 * i] + self[8 * i + 7] >= 0) result++;
        if (partialresult[10] + self[8 * i + 1] + self[8 * i + 2] >= 0) result++;
        if (partialresult[11] + self[8 * i + 1] + self[8 * i + 3] >= 0) result++;
        if (partialresult[12] + self[8 * i + 1] + self[8 * i + 4] >= 0) result++;
        if (partialresult[13] + self[8 * i + 1] + self[8 * i + 5] >= 0) result++;
        if (partialresult[14] + self[8 * i + 1] + self[8 * i + 6] >= 0) result++;
        if (partialresult[15] + self[8 * i + 1] + self[8 * i + 7] >= 0) result++;
        if (partialresult[19] + self[8 * i + 2] + self[8 * i + 3] >= 0) result++;
        if (partialresult[20] + self[8 * i + 2] + self[8 * i + 4] >= 0) result++;
        if (partialresult[21] + self[8 * i + 2] + self[8 * i + 5] >= 0) result++;
        if (partialresult[22] + self[8 * i + 2] + self[8 * i + 6] >= 0) result++;
        if (partialresult[23] + self[8 * i + 2] + self[8 * i + 7] >= 0) result++;
        if (partialresult[28] + self[8 * i + 3] + self[8 * i + 4] >= 0) result++;
        if (partialresult[29] + self[8 * i + 3] + self[8 * i + 5] >= 0) result++;
        if (partialresult[30] + self[8 * i + 3] + self[8 * i + 6] >= 0) result++;
        if (partialresult[31] + self[8 * i + 3] + self[8 * i + 7] >= 0) result++;
        if (partialresult[37] + self[8 * i + 4] + self[8 * i + 5] >= 0) result++;
        if (partialresult[38] + self[8 * i + 4] + self[8 * i + 6] >= 0) result++;
        if (partialresult[39] + self[8 * i + 4] + self[8 * i + 7] >= 0) result++;
        if (partialresult[46] + self[8 * i + 5] + self[8 * i + 6] >= 0) result++;
        if (partialresult[47] + self[8 * i + 5] + self[8 * i + 7] >= 0) result++;
        if (partialresult[55] + self[8 * i + 6] + self[8 * i + 7] >= 0) result++;
    }
    printf("Result after self: %d\n", result);
    // static load balancing
    long long overall_load = 0 ;
    for (int i=0 ; i<n/8 ; i++)
        overall_load += L3[i+nn/8] - L2[i+nn/8] + L1[i+nn/8] - i ;
    int loadstart[NUM_THREADS+2];
    for (int i=0 ; i<=NUM_THREADS ; i++)
        loadstart[i] = 0;
    long long cum_load = 0;
    for (int i=0 ; i<n/8 ; i++){
        cum_load += L3[i+nn/8] - L2[i+nn/8] + L1[i+nn/8] - i ;
        loadstart[cum_load*NUM_THREADS/overall_load+1] = i;
    }
    loadstart[NUM_THREADS] = n/8;
    for (int i=1 ; i<=NUM_THREADS ; i++)
        if (loadstart[i] == 0)
            loadstart[i] = loadstart[i-1] ;
// // output the load sizes
//    for (int i=0 ; i<NUM_THREADS ; i++){
//        cum_load = 0 ;
//        for (int j=loadstart[i]; j<loadstart[i+1] ; j++)
//            cum_load += L3[j+nn/8] - L2[j+nn/8] + L1[j+nn/8] - j ;
//        printf("%d %d %ld\n", i, loadstart[i], cum_load);
//    }
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++) {
        int imin = loadstart[par];//n * par / NUM_THREADS / 8;
        int imax = loadstart[par+1];//n * (par + 1) / NUM_THREADS / 8;
        veci64 resultvec = _mm512_setzero_si512();
        long long localresult = 0;

#ifdef CANOLOOP
        //printf("CANOLOOP %d\n", omp_get_thread_num());
        for (int i = imin; i < imax; i++) for (int j = i + 1; j < L1[i + nn / 8]; j++) {
#else
        int i=0; int j=0; FGF_HILBERT_FOR(i, j, n/8, n/8, i>=imin && i<imax && j>i && j < L1[i+nn/8],
                /*printout_vars(HILLOOP_hilbert, i, j, HILLOOP_level, FURHIL_clevel, HILLOOP_c, FURHIL_lb0, FURHIL_ub0, FURHIL_lb1, FURHIL_ub1)*/
                /*FGF_debug(FURHIL_clevel, FURHIL_lb0, FURHIL_ub0, FURHIL_lb1, FURHIL_ub1) &&*/
                FURHIL_ub0 >= imin && FURHIL_lb0 < imax && FURHIL_ub1-1 > FURHIL_lb0 &&
                FURHIL_lb1 < L1[((FURHIL_ub0-1)>>(FURHIL_clevel+1))+(nn>>(FURHIL_clevel+4))]){
#endif
                register vec vi1 = _mm512_load_pd(array + d * 8 * i);
                register vec vi2 = _mm512_load_pd(array + d * (8 * i + 1));
                register vec vj = _mm512_load_pd(array + d * 8 * j);
                register vec sum1 = vi1 * vj;
                register vec sum9 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1));
                register vec sum2 = vi1 * vj;
                register vec sum10 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2));
                register vec sum3 = vi1 * vj;
                register vec sum11 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3));
                register vec sum4 = vi1 * vj;
                register vec sum12 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4));
                register vec sum5 = vi1 * vj;
                register vec sum13 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5));
                register vec sum6 = vi1 * vj;
                register vec sum14 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6));
                register vec sum7 = vi1 * vj;
                register vec sum15 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7));
                register vec sum8 = vi1 * vj;
                register vec sum16 = vi2 * vj;
                for (register int k = 8; k < d; k += 8) {
                    vi1 = _mm512_load_pd(array + d * 8 * i + k);
                    vi2 = _mm512_load_pd(array + d * (8 * i + 1) + k);
                    vj = _mm512_load_pd(array + d * 8 * j + k);
                    sum1 += vi1 * vj; // _mm512_fmadd_pd
                    sum9 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                    sum2 += vi1 * vj; // _mm512_fmadd_pd
                    sum10 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                    sum3 += vi1 * vj; // _mm512_fmadd_pd
                    sum11 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                    sum4 += vi1 * vj; // _mm512_fmadd_pd
                    sum12 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                    sum5 += vi1 * vj; // _mm512_fmadd_pd
                    sum13 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                    sum6 += vi1 * vj; // _mm512_fmadd_pd
                    sum14 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                    sum7 += vi1 * vj; // _mm512_fmadd_pd
                    sum15 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                    sum8 += vi1 * vj;
                    sum16 += vi2 * vj;
                }
                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
                vec hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                // _mm512_store_pd(partialresult, hhh) ;
                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
                hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 1]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                vi1 = _mm512_load_pd(array + d * (8 * i + 2));
                vi2 = _mm512_load_pd(array + d * (8 * i + 3));
                vj = _mm512_load_pd(array + d * 8 * j);
                sum1 = vi1 * vj;
                sum9 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1));
                sum2 = vi1 * vj;
                sum10 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2));
                sum3 = vi1 * vj;
                sum11 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3));
                sum4 = vi1 * vj;
                sum12 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4));
                sum5 = vi1 * vj;
                sum13 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5));
                sum6 = vi1 * vj;
                sum14 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6));
                sum7 = vi1 * vj;
                sum15 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7));
                sum8 = vi1 * vj;
                sum16 = vi2 * vj;
                for (register int k = 8; k < d; k += 8) {
                    vi1 = _mm512_load_pd(array + d * (8 * i + 2) + k);
                    vi2 = _mm512_load_pd(array + d * (8 * i + 3) + k);
                    vj = _mm512_load_pd(array + d * 8 * j + k);
                    sum1 += vi1 * vj; // _mm512_fmadd_pd
                    sum9 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                    sum2 += vi1 * vj; // _mm512_fmadd_pd
                    sum10 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                    sum3 += vi1 * vj; // _mm512_fmadd_pd
                    sum11 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                    sum4 += vi1 * vj; // _mm512_fmadd_pd
                    sum12 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                    sum5 += vi1 * vj; // _mm512_fmadd_pd
                    sum13 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                    sum6 += vi1 * vj; // _mm512_fmadd_pd
                    sum14 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                    sum7 += vi1 * vj; // _mm512_fmadd_pd
                    sum15 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                    sum8 += vi1 * vj;
                    sum16 += vi2 * vj;
                }
                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
                hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 2]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                // _mm512_store_pd(partialresult, hhh) ;
                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
                hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 3]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                vi1 = _mm512_load_pd(array + d * (8 * i + 4));
                vi2 = _mm512_load_pd(array + d * (8 * i + 5));
                vj = _mm512_load_pd(array + d * 8 * j);
                sum1 = vi1 * vj;
                sum9 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1));
                sum2 = vi1 * vj;
                sum10 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2));
                sum3 = vi1 * vj;
                sum11 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3));
                sum4 = vi1 * vj;
                sum12 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4));
                sum5 = vi1 * vj;
                sum13 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5));
                sum6 = vi1 * vj;
                sum14 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6));
                sum7 = vi1 * vj;
                sum15 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7));
                sum8 = vi1 * vj;
                sum16 = vi2 * vj;
                for (register int k = 8; k < d; k += 8) {
                    vi1 = _mm512_load_pd(array + d * (8 * i + 4) + k);
                    vi2 = _mm512_load_pd(array + d * (8 * i + 5) + k);
                    vj = _mm512_load_pd(array + d * 8 * j + k);
                    sum1 += vi1 * vj; // _mm512_fmadd_pd
                    sum9 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                    sum2 += vi1 * vj; // _mm512_fmadd_pd
                    sum10 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                    sum3 += vi1 * vj; // _mm512_fmadd_pd
                    sum11 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                    sum4 += vi1 * vj; // _mm512_fmadd_pd
                    sum12 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                    sum5 += vi1 * vj; // _mm512_fmadd_pd
                    sum13 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                    sum6 += vi1 * vj; // _mm512_fmadd_pd
                    sum14 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                    sum7 += vi1 * vj; // _mm512_fmadd_pd
                    sum15 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                    sum8 += vi1 * vj;
                    sum16 += vi2 * vj;
                }
                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
                hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 4]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                // _mm512_store_pd(partialresult, hhh) ;
                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
                hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 5]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                vi1 = _mm512_load_pd(array + d * (8 * i + 6));
                vi2 = _mm512_load_pd(array + d * (8 * i + 7));
                vj = _mm512_load_pd(array + d * 8 * j);
                sum1 = vi1 * vj;
                sum9 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1));
                sum2 = vi1 * vj;
                sum10 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2));
                sum3 = vi1 * vj;
                sum11 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3));
                sum4 = vi1 * vj;
                sum12 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4));
                sum5 = vi1 * vj;
                sum13 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5));
                sum6 = vi1 * vj;
                sum14 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6));
                sum7 = vi1 * vj;
                sum15 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7));
                sum8 = vi1 * vj;
                sum16 = vi2 * vj;
                for (register int k = 8; k < d; k += 8) {
                    vi1 = _mm512_load_pd(array + d * (8 * i + 6) + k);
                    vi2 = _mm512_load_pd(array + d * (8 * i + 7) + k);
                    vj = _mm512_load_pd(array + d * 8 * j + k);
                    sum1 += vi1 * vj; // _mm512_fmadd_pd
                    sum9 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                    sum2 += vi1 * vj; // _mm512_fmadd_pd
                    sum10 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                    sum3 += vi1 * vj; // _mm512_fmadd_pd
                    sum11 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                    sum4 += vi1 * vj; // _mm512_fmadd_pd
                    sum12 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                    sum5 += vi1 * vj; // _mm512_fmadd_pd
                    sum13 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                    sum6 += vi1 * vj; // _mm512_fmadd_pd
                    sum14 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                    sum7 += vi1 * vj; // _mm512_fmadd_pd
                    sum15 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                    sum8 += vi1 * vj;
                    sum16 += vi2 * vj;
                }
                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
                hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 6]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                // _mm512_store_pd(partialresult, hhh) ;
                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
                hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 7]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                localresult += 64;
//            }  else {
//                printf("Shit... %d %d %d %d %d\n", FURHIL_clevel, FURHIL_lb0, FURHIL_ub0, FURHIL_lb1, FURHIL_ub1);
            }
#ifndef CANOLOOP
        FGF_HILBERT_END(i,j);

        i=0; j=0; FGF_HILBERT_FOR(i, j, n/8, n/8, i>=imin && i<imax && j>=L2[i+nn/8] && j < L3[i+nn/8],
                /*printout_vars(HILLOOP_hilbert, i, j, HILLOOP_level, FURHIL_clevel, HILLOOP_c, FURHIL_lb0, FURHIL_ub0, FURHIL_lb1, FURHIL_ub1) &&*/
                FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&
                FURHIL_ub1-1 > L2[(FURHIL_lb0>>(FURHIL_clevel+1))+(nn>>(FURHIL_clevel+4))] &&
                FURHIL_lb1 < L3[((FURHIL_ub0-1)>>(FURHIL_clevel+1))+(nn>>(FURHIL_clevel+4))]){
#else
        for (int i = imin; i < imax; i++) for (int j = L2[i + nn / 8]; j < L3[i + nn / 8]; j++) {
#endif
                register vec vi1 = _mm512_load_pd(array + d * 8 * i);
                register vec vi2 = _mm512_load_pd(array + d * (8 * i + 1));
                register vec vj = _mm512_load_pd(array + d * 8 * j);
                register vec sum1 = vi1 * vj;
                register vec sum9 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1));
                register vec sum2 = vi1 * vj;
                register vec sum10 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2));
                register vec sum3 = vi1 * vj;
                register vec sum11 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3));
                register vec sum4 = vi1 * vj;
                register vec sum12 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4));
                register vec sum5 = vi1 * vj;
                register vec sum13 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5));
                register vec sum6 = vi1 * vj;
                register vec sum14 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6));
                register vec sum7 = vi1 * vj;
                register vec sum15 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7));
                register vec sum8 = vi1 * vj;
                register vec sum16 = vi2 * vj;
                for (register int k = 8; k < d; k += 8) {
                    vi1 = _mm512_load_pd(array + d * 8 * i + k);
                    vi2 = _mm512_load_pd(array + d * (8 * i + 1) + k);
                    vj = _mm512_load_pd(array + d * 8 * j + k);
                    sum1 += vi1 * vj; // _mm512_fmadd_pd
                    sum9 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                    sum2 += vi1 * vj; // _mm512_fmadd_pd
                    sum10 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                    sum3 += vi1 * vj; // _mm512_fmadd_pd
                    sum11 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                    sum4 += vi1 * vj; // _mm512_fmadd_pd
                    sum12 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                    sum5 += vi1 * vj; // _mm512_fmadd_pd
                    sum13 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                    sum6 += vi1 * vj; // _mm512_fmadd_pd
                    sum14 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                    sum7 += vi1 * vj; // _mm512_fmadd_pd
                    sum15 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                    sum8 += vi1 * vj;
                    sum16 += vi2 * vj;
                }
                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
                vec hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                // _mm512_store_pd(partialresult, hhh) ;
                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
                hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 1]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                vi1 = _mm512_load_pd(array + d * (8 * i + 2));
                vi2 = _mm512_load_pd(array + d * (8 * i + 3));
                vj = _mm512_load_pd(array + d * 8 * j);
                sum1 = vi1 * vj;
                sum9 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1));
                sum2 = vi1 * vj;
                sum10 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2));
                sum3 = vi1 * vj;
                sum11 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3));
                sum4 = vi1 * vj;
                sum12 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4));
                sum5 = vi1 * vj;
                sum13 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5));
                sum6 = vi1 * vj;
                sum14 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6));
                sum7 = vi1 * vj;
                sum15 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7));
                sum8 = vi1 * vj;
                sum16 = vi2 * vj;
                for (register int k = 8; k < d; k += 8) {
                    vi1 = _mm512_load_pd(array + d * (8 * i + 2) + k);
                    vi2 = _mm512_load_pd(array + d * (8 * i + 3) + k);
                    vj = _mm512_load_pd(array + d * 8 * j + k);
                    sum1 += vi1 * vj; // _mm512_fmadd_pd
                    sum9 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                    sum2 += vi1 * vj; // _mm512_fmadd_pd
                    sum10 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                    sum3 += vi1 * vj; // _mm512_fmadd_pd
                    sum11 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                    sum4 += vi1 * vj; // _mm512_fmadd_pd
                    sum12 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                    sum5 += vi1 * vj; // _mm512_fmadd_pd
                    sum13 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                    sum6 += vi1 * vj; // _mm512_fmadd_pd
                    sum14 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                    sum7 += vi1 * vj; // _mm512_fmadd_pd
                    sum15 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                    sum8 += vi1 * vj;
                    sum16 += vi2 * vj;
                }
                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
                hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 2]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                // _mm512_store_pd(partialresult, hhh) ;
                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
                hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 3]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                vi1 = _mm512_load_pd(array + d * (8 * i + 4));
                vi2 = _mm512_load_pd(array + d * (8 * i + 5));
                vj = _mm512_load_pd(array + d * 8 * j);
                sum1 = vi1 * vj;
                sum9 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1));
                sum2 = vi1 * vj;
                sum10 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2));
                sum3 = vi1 * vj;
                sum11 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3));
                sum4 = vi1 * vj;
                sum12 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4));
                sum5 = vi1 * vj;
                sum13 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5));
                sum6 = vi1 * vj;
                sum14 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6));
                sum7 = vi1 * vj;
                sum15 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7));
                sum8 = vi1 * vj;
                sum16 = vi2 * vj;
                for (register int k = 8; k < d; k += 8) {
                    vi1 = _mm512_load_pd(array + d * (8 * i + 4) + k);
                    vi2 = _mm512_load_pd(array + d * (8 * i + 5) + k);
                    vj = _mm512_load_pd(array + d * 8 * j + k);
                    sum1 += vi1 * vj; // _mm512_fmadd_pd
                    sum9 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                    sum2 += vi1 * vj; // _mm512_fmadd_pd
                    sum10 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                    sum3 += vi1 * vj; // _mm512_fmadd_pd
                    sum11 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                    sum4 += vi1 * vj; // _mm512_fmadd_pd
                    sum12 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                    sum5 += vi1 * vj; // _mm512_fmadd_pd
                    sum13 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                    sum6 += vi1 * vj; // _mm512_fmadd_pd
                    sum14 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                    sum7 += vi1 * vj; // _mm512_fmadd_pd
                    sum15 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                    sum8 += vi1 * vj;
                    sum16 += vi2 * vj;
                }
                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
                hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 4]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                // _mm512_store_pd(partialresult, hhh) ;
                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
                hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 5]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                vi1 = _mm512_load_pd(array + d * (8 * i + 6));
                vi2 = _mm512_load_pd(array + d * (8 * i + 7));
                vj = _mm512_load_pd(array + d * 8 * j);
                sum1 = vi1 * vj;
                sum9 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 1));
                sum2 = vi1 * vj;
                sum10 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 2));
                sum3 = vi1 * vj;
                sum11 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 3));
                sum4 = vi1 * vj;
                sum12 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 4));
                sum5 = vi1 * vj;
                sum13 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 5));
                sum6 = vi1 * vj;
                sum14 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 6));
                sum7 = vi1 * vj;
                sum15 = vi2 * vj;
                vj = _mm512_load_pd(array + d * (8 * j + 7));
                sum8 = vi1 * vj;
                sum16 = vi2 * vj;
                for (register int k = 8; k < d; k += 8) {
                    vi1 = _mm512_load_pd(array + d * (8 * i + 6) + k);
                    vi2 = _mm512_load_pd(array + d * (8 * i + 7) + k);
                    vj = _mm512_load_pd(array + d * 8 * j + k);
                    sum1 += vi1 * vj; // _mm512_fmadd_pd
                    sum9 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                    sum2 += vi1 * vj; // _mm512_fmadd_pd
                    sum10 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                    sum3 += vi1 * vj; // _mm512_fmadd_pd
                    sum11 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                    sum4 += vi1 * vj; // _mm512_fmadd_pd
                    sum12 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                    sum5 += vi1 * vj; // _mm512_fmadd_pd
                    sum13 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                    sum6 += vi1 * vj; // _mm512_fmadd_pd
                    sum14 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                    sum7 += vi1 * vj; // _mm512_fmadd_pd
                    sum15 += vi2 * vj;
                    vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                    sum8 += vi1 * vj;
                    sum16 += vi2 * vj;
                }
                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
                hhh = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 6]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                // _mm512_store_pd(partialresult, hhh) ;
                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
                hhh = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[8 * i + 7]);
                resultvec += _mm512_srli_epi64(_mm512_castpd_si512(hhh), 63);
                localresult += 64;
            }
#ifndef CANOLOOP
        FGF_HILBERT_END(i,j);
#endif
        long long resultvecX[8];
        _mm512_storeu_si512((void *) resultvecX, resultvec);
        //    printf("%ld %ld %ld %ld %ld %ld %ld %ld\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
        //    printf("%lx %lx %lx %lx %lx %lx %lx %lx\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
        localresult -= resultvecX[0] + resultvecX[1] + resultvecX[2] + resultvecX[3] + resultvecX[4] + resultvecX[5] + resultvecX[6] + resultvecX[7];
#pragma omp critical
        {
            result += localresult;
        }
    }
    printf("Results nlj: %d\n", 2 * result + n);
    return result;
}

int joinCanonic(int n, int d, double epsilon, double* array) {
    int result = 0;
    double epsilon2 = epsilon*epsilon;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dist = 0.;
            for (int k = 0; k < d; k++) {
                double h = array[i * d + k] - array[j * d + k];
                dist += h*h;
            }
            if (dist <= epsilon2)
                result++;
        }
    }
    printf("Results nlj: %d\n", 2 * result + n);
    printf("Join nlj   | %8.3f | %8.3f\n", stopc(), stop());
    return 2 * result + n;
}

int joinCanonicVerbose(int n, int d, double epsilon, double* array) {
    int result = 0;
    double epsilon2 = epsilon*epsilon;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double dist = 0.;
            for (int k = 0; k < d; k++) {
                double h = array[i * d + k] - array[j * d + k];
                dist += h*h;
            }
            printf("%4.2f", dist);
            if (dist <= epsilon2) {
                result++;
                printf("* ");
            } else printf("  ");
        }
        printf("\n");
    }
    printf("\n");
    return 2 * result + n;
}

void prepareTwoStripes(int n, int d, double epsilon, double *array, int ** lower, int ** upper, double *self){
    epsilonGridOrdering(n, d, epsilon, array);
    int nn = ceilpowtwo(n);
    lower[0] = (int *) mallocA64(sizeof (int) * 2 * nn);
    upper[0] = (int *) mallocA64(sizeof (int) * 2 * nn);
    lower[1] = (int *) mallocA64(sizeof (int) * 2 * nn);
    upper[1] = (int *) mallocA64(sizeof (int) * 2 * nn);
    for (int i=0 ; i<n ; i++)
        lower[0][i+nn] = i ;
    epsilonGridFillList3(n, d, epsilon, array, upper[0] + nn, lower[1] + nn, upper[1] + nn);
//    for (int i=0 ; i<20 ; i++)
//        printf("i=%d >>> %d %d %d\n",i, upper[0][nn+i], lower[1][nn+i], upper[1][nn+i]);
    for (int j=0 ; j<2 ; j++){
        for (int i = 0; i < n; i++)
           lower[j][nn + i] = lower[j][nn + i] / 8;
        for (int i = 0; i < n; i++)
            upper[j][nn + i] = (upper[j][nn + i] + 7) / 8;
        epsilonGridCompleteListMin(nn, lower[j]);
        epsilonGridCompleteListMax(nn, upper[j]);
    }
    double epsilon22 = epsilon*epsilon/2;
    for (int i=0 ; i<n ; i++){
        double h=epsilon22;
        for(int j=0 ; j<d ; j++)
            h-=array[i*d+j]*array[i*d+j];
        self[i]=h/2;
    }
}
//char testarray[25000][25000];

int epsilonGridJoinSimplified(int n, int d, double epsilon, double* array, int stripes, int ** lower, int ** upper, double *self) {
    double value = -1.;
    long long result = 0;
    double epsilon2 = epsilon * epsilon;
    veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};
    veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};
    veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};
    veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};
    veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};
    veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};

    long long overall_load = 0;
    for (int i = 0; i < n / 8; i++)
        for (int j = 0; j < stripes; j++)
            overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];
    int loadstart[NUM_THREADS + 2];
    for (int i = 0; i <= NUM_THREADS; i++)
        loadstart[i] = 0;
    long long cum_load = 0;
    for (int i = 0; i < n / 8; i++) {
        for (int j = 0; j < stripes; j++)
            cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];
        loadstart[cum_load * NUM_THREADS / overall_load + 1] = i;
    }
    loadstart[NUM_THREADS] = n / 8;
    for (int i = 1; i <= NUM_THREADS; i++)
        if (loadstart[i] == 0)
            loadstart[i] = loadstart[i - 1];
    // // output the load sizes
    //    for (int i=0 ; i<NUM_THREADS ; i++){
    //        cum_load = 0 ;
    //        for (int j=loadstart[i]; j<loadstart[i+1] ; j++)
    //            cum_load += L3[j+nn/8] - L2[j+nn/8] + L1[j+nn/8] - j ;
    //        printf("%d %d %ld\n", i, loadstart[i], cum_load);
    //    }
#pragma omp parallel for reduction(+:result)
    for (int par = 0; par < NUM_THREADS; par++) {
        int imin = loadstart[par]; //n * par / NUM_THREADS / 8;
        int imax = loadstart[par + 1]; //n * (par + 1) / NUM_THREADS / 8;
        veci64 resultvec = _mm512_setzero_si512();
//        long long localresult = 0;

        for (int s = 0; s < stripes; s++) {
#ifndef CANOLOOP
            int i = 0;
            int j = 0;

            FGF_HILBERT_FOR(i, j, n / 8, n / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],
                    /*printout_vars(HILLOOP_hilbert, i, j, HILLOOP_level, FURHIL_clevel, HILLOOP_c, FURHIL_lb0, FURHIL_ub0, FURHIL_lb1, FURHIL_ub1) &&*/
                    FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&
                    FURHIL_ub1 - 1 >= lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&
                    FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {
//            FGF_HILBERT_FOR(i, j, n / 8, n / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],
//                    /*printout_vars(HILLOOP_hilbert, i, j, HILLOOP_level, FURHIL_clevel, HILLOOP_c, FURHIL_lb0, FURHIL_ub0, FURHIL_lb1, FURHIL_ub1) &&*/
//                    FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&
//                    FURHIL_ub1 - 1 > lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&
//                    FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {
#else
            for (int i = imin; i < imax; i++) for (int j = L2[i + nn / 8]; j < L3[i + nn / 8]; j++) {
#endif
//                testarray[i][j] ++ ;
                    for (int II = 8 * i; II < 8 * i + 8; II += 2) {
                        register vec vi1 = _mm512_load_pd(array + d * II);
                        register vec vi2 = _mm512_load_pd(array + d * (II + 1));
                        register vec vj = _mm512_load_pd(array + d * 8 * j);
                        register vec sum1 = _mm512_mul_pd(vi1, vj);
                        register vec sum9 = _mm512_mul_pd(vi2, vj);
                        vj = _mm512_load_pd(array + d * (8 * j + 1));
                        register vec sum2 = _mm512_mul_pd(vi1, vj);
                        register vec sum10 = _mm512_mul_pd(vi2, vj);
                        vj = _mm512_load_pd(array + d * (8 * j + 2));
                        register vec sum3 = _mm512_mul_pd(vi1, vj);
                        register vec sum11 = _mm512_mul_pd(vi2, vj);
                        vj = _mm512_load_pd(array + d * (8 * j + 3));
                        register vec sum4 = _mm512_mul_pd(vi1, vj);
                        register vec sum12 = _mm512_mul_pd(vi2, vj);
                        vj = _mm512_load_pd(array + d * (8 * j + 4));
                        register vec sum5 = _mm512_mul_pd(vi1, vj);
                        register vec sum13 = _mm512_mul_pd(vi2, vj);
                        vj = _mm512_load_pd(array + d * (8 * j + 5));
                        register vec sum6 = _mm512_mul_pd(vi1, vj);
                        register vec sum14 = _mm512_mul_pd(vi2, vj);
                        vj = _mm512_load_pd(array + d * (8 * j + 6));
                        register vec sum7 = _mm512_mul_pd(vi1, vj);
                        register vec sum15 = _mm512_mul_pd(vi2, vj);
                        vj = _mm512_load_pd(array + d * (8 * j + 7));
                        register vec sum8 = _mm512_mul_pd(vi1, vj);
                        register vec sum16 = _mm512_mul_pd(vi2, vj);
                        for (register int k = 8; k < d; k += 8) {
                            vi1 = _mm512_load_pd(array + d * II + k);
                            vi2 = _mm512_load_pd(array + d * (II + 1) + k);
                            vj = _mm512_load_pd(array + d * 8 * j + k);
                            sum1 = _mm512_fmadd_pd(vi1, vj, sum1);
                            sum9 = _mm512_fmadd_pd(vi2, vj, sum9);
                            vj = _mm512_load_pd(array + d * (8 * j + 1) + k);
                            sum2 = _mm512_fmadd_pd(vi1, vj, sum2); // _mm512_fmadd_pd
                            sum10 = _mm512_fmadd_pd(vi2, vj, sum10);
                            vj = _mm512_load_pd(array + d * (8 * j + 2) + k);
                            sum3 = _mm512_fmadd_pd(vi1, vj, sum3); // _mm512_fmadd_pd
                            sum11 = _mm512_fmadd_pd(vi2, vj, sum11);
                            vj = _mm512_load_pd(array + d * (8 * j + 3) + k);
                            sum4 = _mm512_fmadd_pd(vi1, vj, sum4); // _mm512_fmadd_pd
                            sum12 = _mm512_fmadd_pd(vi2, vj, sum12);
                            vj = _mm512_load_pd(array + d * (8 * j + 4) + k);
                            sum5 = _mm512_fmadd_pd(vi1, vj, sum5); // _mm512_fmadd_pd
                            sum13 = _mm512_fmadd_pd(vi2, vj, sum13);
                            vj = _mm512_load_pd(array + d * (8 * j + 5) + k);
                            sum6 = _mm512_fmadd_pd(vi1, vj, sum6); // _mm512_fmadd_pd
                            sum14 = _mm512_fmadd_pd(vi2, vj, sum14);
                            vj = _mm512_load_pd(array + d * (8 * j + 6) + k);
                            sum7 = _mm512_fmadd_pd(vi1, vj, sum7); // _mm512_fmadd_pd
                            sum15 = _mm512_fmadd_pd(vi2, vj, sum15);
                            vj = _mm512_load_pd(array + d * (8 * j + 7) + k);
                            sum8 = _mm512_fmadd_pd(vi1, vj, sum8);
                            sum16 = _mm512_fmadd_pd(vi2, vj, sum16);
                        }
                        for (int III = 0; III < 2; III++) {
                            vec indicator;
                            if (III) {
                                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);
                                indicator = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +
                                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[II + 1]);
                                if(i==j)
                                    indicator = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(indicator), __mmask8(255 >> (7-(II&7|1))), *(long long *)&value));
                            } else {
                                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);
                                indicator = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +
                                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[II]);
                                if(i==j)
                                    indicator = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(indicator), __mmask8(255 >> (7-(II&6))), *(long long *)&value));
                            }
                            // THIS IS THE ACTUAL LOOP BODY OF THE SIMILARITY JOIN
                            resultvec += _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(indicator), 63));
//                            resultvec += _mm512_srli_epi64(_mm512_castpd_si512(indicator), 63);
                            // UNTIL HERE (THIS IS THE ACTUAL LOOP BODY OF THE SIMILARITY JOIN)
                        }
                    }

#ifndef CANOLOOP
                }
                FGF_HILBERT_END(i, j);
#else
                }
#endif
        }
        long long resultvecX[8];
        _mm512_storeu_si512((void *) resultvecX, resultvec);
        //    printf("%ld %ld %ld %ld %ld %ld %ld %ld\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
        //    printf("%lx %lx %lx %lx %lx %lx %lx %lx\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
        result += resultvecX[0] + resultvecX[1] + resultvecX[2] + resultvecX[3] + resultvecX[4] + resultvecX[5] + resultvecX[6] + resultvecX[7];
//#pragma omp critical
//        {
//            result += localresult;
//        }
    }
    printf("Results: %d\n", 2 * result + n);
    return result;
}

# define EGO_PARALLEL(n,d,epsilon,array) {\
    int EGO_n = (n);\
    int EGO_d = (d);\
    double EGO_epsilon = (epsilon);\
    double * EGO_array = (array);\
    int stripes = 2;\
    epsilonGridOrdering(EGO_n, EGO_d, EGO_epsilon, EGO_array);\
    int nn = ceilpowtwo(EGO_n);\
    int **lower = (int **) malloc (2*sizeof(int*));\
    int **upper = (int **) malloc (2*sizeof(int*));\
    double *self = mallocA64(sizeof (double) * EGO_n);\
    lower[0] = (int *) mallocA64(sizeof (int) * 2 * nn);\
    upper[0] = (int *) mallocA64(sizeof (int) * 2 * nn);\
    lower[1] = (int *) mallocA64(sizeof (int) * 2 * nn);\
    upper[1] = (int *) mallocA64(sizeof (int) * 2 * nn);\
    for (int i=0 ; i<EGO_n ; i++)\
        lower[0][i+nn] = i ;\
    epsilonGridFillList3(EGO_n, EGO_d, EGO_epsilon, EGO_array, upper[0] + nn, lower[1] + nn, upper[1] + nn);\
    for (int j=0 ; j<2 ; j++){\
        for (int i = 0; i < EGO_n; i++)\
           lower[j][nn + i] = lower[j][nn + i] / 8;\
        for (int i = 0; i < EGO_n; i++)\
            upper[j][nn + i] = (upper[j][nn + i] + 7) / 8;\
        epsilonGridCompleteListMin(nn, lower[j]);\
        epsilonGridCompleteListMax(nn, upper[j]);\
    }\
    EGO_epsilon = EGO_epsilon * EGO_epsilon / 2;\
    for (int i=0 ; i<EGO_n ; i++){\
        double h=EGO_epsilon;\
        for(int j=0 ; j<EGO_d ; j++)\
            h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];\
        self[i]=h/2;\
    }\
    double value = -1.;\
    veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};\
    veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};\
    veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};\
    veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};\
    veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};\
    veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};\
    long long overall_load = 0;\
    for (int i = 0; i < EGO_n / 8; i++)\
        for (int j = 0; j < stripes; j++)\
            overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
    int loadstart[NUM_THREADS + 2];\
    for (int i = 0; i <= NUM_THREADS; i++)\
        loadstart[i] = 0;\
    long long cum_load = 0;\
    for (int i = 0; i < EGO_n / 8; i++) {\
        for (int j = 0; j < stripes; j++)\
            cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
        loadstart[cum_load * NUM_THREADS / overall_load + 1] = i;\
    }\
    loadstart[NUM_THREADS] = EGO_n / 8;\
    for (int i = 1; i <= NUM_THREADS; i++)\
        if (loadstart[i] == 0)\
            loadstart[i] = loadstart[i - 1];

#define EGO_PREPARE for (int par = 0; par < NUM_THREADS; par++) {

#define EGO_LOOP\
        int imin = loadstart[par];\
        int imax = loadstart[par + 1];\
        for (int s = 0; s < stripes; s++) {\
            int i = 0;\
            int j = 0;\
            FGF_HILBERT_FOR(i, j, EGO_n / 8, EGO_n / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],\
                    FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&\
                    FURHIL_ub1 - 1 > lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&\
                    FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {\
                    for (int II = 8 * i; II < 8 * i + 8; II += 2) {\
                        register vec vi1 = _mm512_load_pd(EGO_array + EGO_d * II);\
                        register vec vi2 = _mm512_load_pd(EGO_array + EGO_d * (II + 1));\
                        register vec vj = _mm512_load_pd(EGO_array + EGO_d * 8 * j);\
                        register vec sum1 = _mm512_mul_pd(vi1, vj);\
                        register vec sum9 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 1));\
                        register vec sum2 = _mm512_mul_pd(vi1, vj);\
                        register vec sum10 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 2));\
                        register vec sum3 = _mm512_mul_pd(vi1, vj);\
                        register vec sum11 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 3));\
                        register vec sum4 = _mm512_mul_pd(vi1, vj);\
                        register vec sum12 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 4));\
                        register vec sum5 = _mm512_mul_pd(vi1, vj);\
                        register vec sum13 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 5));\
                        register vec sum6 = _mm512_mul_pd(vi1, vj);\
                        register vec sum14 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 6));\
                        register vec sum7 = _mm512_mul_pd(vi1, vj);\
                        register vec sum15 = _mm512_mul_pd(vi2, vj);\
                        vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 7));\
                        register vec sum8 = _mm512_mul_pd(vi1, vj);\
                        register vec sum16 = _mm512_mul_pd(vi2, vj);\
                        for (register int k = 8; k < d; k += 8) {\
                            vi1 = _mm512_load_pd(EGO_array + EGO_d * II + k);\
                            vi2 = _mm512_load_pd(EGO_array + EGO_d * (II + 1) + k);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * 8 * j + k);\
                            sum1 = _mm512_fmadd_pd(vi1, vj, sum1); \
                            sum9 = _mm512_fmadd_pd(vi2, vj, sum9);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 1) + k);\
                            sum2 = _mm512_fmadd_pd(vi1, vj, sum2);\
                            sum10 = _mm512_fmadd_pd(vi2, vj, sum10);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 2) + k);\
                            sum3 = _mm512_fmadd_pd(vi1, vj, sum3);\
                            sum11 = _mm512_fmadd_pd(vi2, vj, sum11);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 3) + k);\
                            sum4 = _mm512_fmadd_pd(vi1, vj, sum4);\
                            sum12 = _mm512_fmadd_pd(vi2, vj, sum12);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 4) + k);\
                            sum5 = _mm512_fmadd_pd(vi1, vj, sum5);\
                            sum13 = _mm512_fmadd_pd(vi2, vj, sum13);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 5) + k);\
                            sum6 = _mm512_fmadd_pd(vi1, vj, sum6);\
                            sum14 = _mm512_fmadd_pd(vi2, vj, sum14);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 6) + k);\
                            sum7 = _mm512_fmadd_pd(vi1, vj, sum7);\
                            sum15 = _mm512_fmadd_pd(vi2, vj, sum15);\
                            vj = _mm512_load_pd(EGO_array + EGO_d * (8 * j + 7) + k);\
                            sum8 = _mm512_fmadd_pd(vi1, vj, sum8);\
                            sum16 = _mm512_fmadd_pd(vi2, vj, sum16);\
                        }\
                        for (int III = 0; III < 2; III++) {\
                            vec indicator;\
                            if (III) {\
                                transposeAVX512(sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16);\
                                indicator = sum9 + sum10 + sum11 + sum12 + sum13 + sum14 + sum15 + sum16 +\
                                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[II + 1]);\
                                if(i==j)\
                                    indicator = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(indicator), __mmask8(255 >> (7-(II&7|1))), *(long long *)&value));\
                            } else {\
                                transposeAVX512(sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8);\
                                indicator = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8 +\
                                        _mm512_load_pd(self + j * 8) + _mm512_set1_pd(self[II]);\
                                if(i==j)\
                                    indicator = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(indicator), __mmask8(255 >> (7-(II&6))), *(long long *)&value));\
                            }

#ifndef KBLOCK
#define KBLOCK 16
#endif

void prepareStripes(int n, int d, int numStripes, double epsilon, double *array, int ** lower, int ** upper, double *self){
//                    for (int a=0 ; a<d && a<20 ; a++)
//                        printf("%5.3lf ", array[17*1024*d+a]);
//                    printf("\n");
//                    for (int a=0 ; a<d && a<d ; a++)
//                        printf("%3.1lf ", (double)(int)(array[17*1024*d+a]/epsilon)*epsilon);
//                    printf("\n");
    double lowerDisplace[numStripes][d];
    double upperDisplace[numStripes][d];
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
        printf("%d %d %d\n", t, i, umrechnung);
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

void prepareStripesFull(int n, int d, int numStripes, double epsilon, double *array, int ** lower, int ** upper){
    double lowerDisplace[numStripes][d];
    double upperDisplace[numStripes][d];
    int neun = 1;
    int activeDimensions = 0;
    while (neun < numStripes){
        neun*=3 ;
        activeDimensions++;
    }
    printf("activeDimensions: %d, %d\n", activeDimensions, neun);
    for (int t=0 ; t<numStripes-1 ; t++){
        int k=t+1;
        int i=k*neun/numStripes;
        int i2=(k-1)*neun/numStripes;
        if ((i + 3) / 3 > (i2 + 3) / 3)
            i -= i % 3;
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
        lowerDisplace[0][j] = (-epsilon);
    for(int j=0; j<d ; j++)
        upperDisplace[numStripes-1][j] = epsilon;
    for(int i=0 ; i<numStripes ; i++){
        printf("%5.2lf %5.2lf %5.2lf %5.2lf %5.2lf     %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf\n",
                lowerDisplace[i][0],lowerDisplace[i][1],lowerDisplace[i][2],lowerDisplace[i][3],lowerDisplace[i][4],
                upperDisplace[i][0],upperDisplace[i][1],upperDisplace[i][2],upperDisplace[i][3],upperDisplace[i][4]) ;
    }
    lower[0] = (int *) callocA64(sizeof (int) * 2 * n * numStripes);
    for (int j=0 ; j<numStripes ; j++){
        lower[j] = (int*)((char *)(lower[0]) + j  * 2 * n * sizeof(int)) ;
        upper[j] = (int*)((char *)(lower[j]) +  n * sizeof(int)) ;
    }
#pragma omp parallel for
    for (int par = 0; par < NUM_THREADS; par++){
        int imin = par * n / NUM_THREADS / 8 * 8;
        int imax = (par + 1) * n / NUM_THREADS / 8 * 8;
        if (par+1 == NUM_THREADS)
            imax = n;
        double h[d];
        for(int j=0; j<numStripes ; j++){
            for(int a=0 ; a<d ; a++)
                h[a] = array[imin*d+a]+lowerDisplace[j][a] ;
            int a = 0 ;
            int b = n-1 ;
            int m = (a+b)/2;
            while (b-a > 1){
                if(epsilonGridCompare(array + m*d, h) >= 0)
                    b = m ;
                else
                    a = m ;
                m = (a+b)/2 ;
            }
            for(int i=imin ; i<imax ; i++) {
                for(int a=0 ; a<d ; a++)
                    h[a] = array[i*d+a]+lowerDisplace[j][a] ;
                while (m > 0 && epsilonGridCompare(array + m * d, h) >= 0)
                    m--;
                while (m < n && epsilonGridCompare(array + m * d, h) < 0)
                    m++;
                lower[j][i] = m ;
            }
        }
        for(int j=0; j<numStripes ; j++){
            for(int a=0 ; a<d ; a++)
                h[a] = array[imin*d+a]+upperDisplace[j][a] ;
            int a = 0 ;
            int b = n-1 ;
            int m = (a+b)/2;
            while (b-a > 1){
                if(epsilonGridCompare(array + m*d, h) >= 0)
                    b = m ;
                else
                    a = m ;
                m = (a+b)/2 ;
            }
            for(int i=imin ; i<imax ; i++) {
                for(int a=0 ; a<d ; a++)
                    h[a] = array[i*d+a]+upperDisplace[j][a] ;
                while (m > 0 && epsilonGridCompare(array + m * d, h) >= 0)
                    m--;
                while (m < n && epsilonGridCompare(array + m * d, h) < 0)
                    m++;
                upper[j][i] = m ;
            }
        }
        for (int j=1 ; j<numStripes ; j++)
            for (int i=imin/8 ; i<imax/8 ; i++)
                upper[j-1][i] = min(upper[j-1][i], lower[j][i]);
    }
    for(int i=0; i<n; i+=65536){
        for(int j=0; j<numStripes; j++)
            printf("%d %d ", lower[j][i], upper[j][i]);
        printf("\n");
    }
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

static inline void transpose_8xd(int n, int d, double *EGO_array) {
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

static inline void transpose_dx8(int n, int d, double *EGO_array) {
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

# define EGO_PARALLEL_TRAN(n,d,epsilon,stripes,array) {\
    int EGO_n = (n);\
    int EGO_d = (d);\
    double EGO_epsilon = (epsilon);\
    double * EGO_array = (array);\
    int EGO_blocks = (EGO_d + KBLOCK - 1) / KBLOCK;\
    int EGO_stripes = (stripes);\
    epsilonGridOrdering(EGO_n, EGO_d, EGO_epsilon, EGO_array);\
    int nn = ceilpowtwo(EGO_n);\
    int **lower = (int **) malloc (EGO_stripes*sizeof(int*));\
    int **upper = (int **) malloc (EGO_stripes*sizeof(int*));\
    double *self = callocA64(sizeof (double) * EGO_n * EGO_blocks);\
    prepareStripes(EGO_n, EGO_d, EGO_stripes, EGO_epsilon, EGO_array, lower, upper, (double *)0);\
    EGO_epsilon = EGO_epsilon * EGO_epsilon / 2;\
    for (int i=0 ; i<EGO_n ; i++){\
        double h=EGO_epsilon;\
        for(int j=0 ; j<2*KBLOCK && j<EGO_d ; j++)\
            h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];\
        self[i/8*8*EGO_blocks+i%8]=h/2;\
    }\
    for(int k=2 ; k<EGO_blocks ; k++)\
        for (int i=0 ; i<EGO_n ; i++){\
            double h=0;\
            for(int j=k * KBLOCK ; j<(k+1) * KBLOCK && j<EGO_d ; j++)\
            h-=EGO_array[i*EGO_d+j]*EGO_array[i*EGO_d+j];\
        self[(i/8*EGO_blocks+k)*8+i%8]=h/2;\
    }\
    transpose_8xd(EGO_n, EGO_d, EGO_array);\
    printf("end transpose8xd\n");\
    long long overall_load = 0;\
    for (int i = 0; i < EGO_n / 8; i++)\
        for (int j = 0; j < EGO_stripes; j++)\
            overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
    int loadstart[NUM_THREADS + 2];\
    for (int i = 0; i <= NUM_THREADS; i++)\
        loadstart[i] = 0;\
    long long cum_load = 0;\
    for (int i = 0; i < EGO_n / 8; i++) {\
        for (int j = 0; j < EGO_stripes; j++)\
            cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
        loadstart[cum_load * NUM_THREADS / overall_load + 1] = i;\
    }\
    loadstart[NUM_THREADS] = EGO_n / 8;\
    for (int i = 1; i <= NUM_THREADS; i++)\
        if (loadstart[i] == 0)\
            loadstart[i] = loadstart[i - 1];

//    lower[0] = (int *) mallocA64(sizeof (int) * 2 * nn);\
//    upper[0] = (int *) mallocA64(sizeof (int) * 2 * nn);\
//    lower[1] = (int *) mallocA64(sizeof (int) * 2 * nn);\
//    upper[1] = (int *) mallocA64(sizeof (int) * 2 * nn);\
//    for (int i=0 ; i<EGO_n ; i++)\
//        lower[0][i+nn] = i ;\
//    epsilonGridFillList3(EGO_n, EGO_d, EGO_epsilon, EGO_array, upper[0] + nn, lower[1] + nn, upper[1] + nn);\
//    for (int j=0 ; j<2 ; j++){\
//        for (int i = 0; i < EGO_n; i++)\
//           lower[j][nn + i] = lower[j][nn + i] / 8;\
//        for (int i = 0; i < EGO_n; i++)\
//            upper[j][nn + i] = (upper[j][nn + i] + 7) / 8;\
//        epsilonGridCompleteListMin(nn, lower[j]);\
//        epsilonGridCompleteListMax(nn, upper[j]);\
//    }\

inline void print_vector(vec v){
    double h[8]__attribute__((aligned(64)));
    _mm512_store_pd(h, v);
    printf("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n", h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
}
inline void print_matrix(int i, int j, vec v1,vec v2,vec v3,vec v4,vec v5,vec v6,vec v7,vec v8){
    printf("%d %d\n", i, j);
    print_vector(v1);print_vector(v2);print_vector(v3);print_vector(v4);
    print_vector(v5);print_vector(v6);print_vector(v7);print_vector(v8);
    printf("\n");
}



#define EGO_LOOP_TRAN\
        int imin = loadstart[par];\
        int imax = loadstart[par + 1];\
        for (int s = 0; s < EGO_stripes; s++) {\
            int i = 0;\
            int j = 0;\
            register veci64 const1 = _mm512_set1_epi64(1LL) ;\
            register veci64 const2 = _mm512_add_epi64(const1, const1) ;\
            register veci64 const3 = _mm512_add_epi64(const1, const2) ;\
            register veci64 const4 = _mm512_add_epi64(const1, const3) ;\
            register veci64 const5 = _mm512_add_epi64(const1, const4) ;\
            register veci64 const6 = _mm512_add_epi64(const1, const5) ;\
            register veci64 const7 = _mm512_add_epi64(const1, const6) ;\
            register veci64 const0 = _mm512_setzero_si512();\
            FGF_HILBERT_FOR(i, j, EGO_n / 8, EGO_n / 8, i >= imin && i < imax && j >= lower[s][i + nn / 8] && j < upper[s][i + nn / 8],\
                    FURHIL_ub0 >= imin && FURHIL_lb0 < imax &&\
                    FURHIL_ub1 - 1 >= lower[s][(FURHIL_lb0 >> (FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))] &&\
                    FURHIL_lb1 < upper[s][((FURHIL_ub0 - 1)>>(FURHIL_clevel + 1))+(nn >> (FURHIL_clevel + 4))]) {\
                register vec vi = _mm512_load_pd(self + i * 8 * EGO_blocks);\
                register vec vj = _mm512_load_pd(self + j * 8 * EGO_blocks);\
                register vec sum1 = vi + _mm512_permutexvar_pd(const0, vj);\
                register vec sum2 = vi + _mm512_permutexvar_pd(const1, vj);\
                register vec sum3 = vi + _mm512_permutexvar_pd(const2, vj);\
                register vec sum4 = vi + _mm512_permutexvar_pd(const3, vj);\
                register vec sum5 = vi + _mm512_permutexvar_pd(const4, vj);\
                register vec sum6 = vi + _mm512_permutexvar_pd(const5, vj);\
                register vec sum7 = vi + _mm512_permutexvar_pd(const6, vj);\
                register vec sum8 = vi + _mm512_permutexvar_pd(const7, vj);\
                vi = _mm512_load_pd(EGO_array + (EGO_d * i) * 8);\
                vj = _mm512_load_pd(EGO_array + (EGO_d * j) * 8);\
                sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);\
                sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);\
                sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);\
                sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);\
                sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);\
                sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);\
                sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);\
                sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);\
                int k; for (k=1 ; k<d ; k++){\
                    if(k % KBLOCK == 0 && k > KBLOCK){\
                        register veci64 allind = _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63);\
                        allind += _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63);\
                        allind += _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63);\
                        allind += _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63);\
                        allind += _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63);\
                        allind += _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63);\
                        allind += _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63);\
                        allind += _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63);\
                        if(_mm512_reduce_add_epi64(allind) >= 64) {k=d+1; break;}\
                        vi = _mm512_load_pd(self + (i * EGO_blocks + k/KBLOCK) * 8);\
                        vj = _mm512_load_pd(self + (j * EGO_blocks + k/KBLOCK) * 8);\
                        sum1 += vi + _mm512_permutexvar_pd(const0, vj);\
                        sum2 += vi + _mm512_permutexvar_pd(const1, vj);\
                        sum3 += vi + _mm512_permutexvar_pd(const2, vj);\
                        sum4 += vi + _mm512_permutexvar_pd(const3, vj);\
                        sum5 += vi + _mm512_permutexvar_pd(const4, vj);\
                        sum6 += vi + _mm512_permutexvar_pd(const5, vj);\
                        sum7 += vi + _mm512_permutexvar_pd(const6, vj);\
                        sum8 += vi + _mm512_permutexvar_pd(const7, vj);\
                    }\
                    vi = _mm512_load_pd(EGO_array + (EGO_d * i + k) * 8);\
                    vj = _mm512_load_pd(EGO_array + (EGO_d * j + k) * 8);\
                    sum1 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const0, vj), sum1);\
                    sum2 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const1, vj), sum2);\
                    sum3 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const2, vj), sum3);\
                    sum4 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const3, vj), sum4);\
                    sum5 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const4, vj), sum5);\
                    sum6 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const5, vj), sum6);\
                    sum7 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const6, vj), sum7);\
                    sum8 = _mm512_fmadd_pd(vi, _mm512_permutexvar_pd(const7, vj), sum8);\
                }\
                if(k<=d){{\
                    if (i==j){\
                        sum1 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum1), 255, 0xbff0000000000000ull));\
                        sum2 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum2), 254, 0xbff0000000000000ull));\
                        sum3 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum3), 252, 0xbff0000000000000ull));\
                        sum4 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum4), 248, 0xbff0000000000000ull));\
                        sum5 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum5), 240, 0xbff0000000000000ull));\
                        sum6 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum6), 224, 0xbff0000000000000ull));\
                        sum7 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum7), 192, 0xbff0000000000000ull));\
                        sum8 = _mm512_castsi512_pd(_mm512_mask_set1_epi64(_mm512_castpd_si512(sum8), 128, 0xbff0000000000000ull));\
                    }

#define EGO_CONSOLIDATE }}} FGF_HILBERT_END(i, j); }

#define EGO_END    }}
#define EGO_END_TRAN    }transpose_dx8(EGO_n, EGO_d, EGO_array);}
#define EGO_END_RESUME  }
#define EGO_PARALLEL_TRAN_RESUME

int test_ego_loop(int n, int d, double epsilon, double *array){
    long long result = 0;
    EGO_PARALLEL(n, d, epsilon, array)
        #pragma omp parallel for reduction(+:result)
    EGO_PREPARE
        veci64 resultvec = _mm512_setzero_si512();
    EGO_LOOP{
        resultvec += _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(indicator), 63));
    }
    EGO_CONSOLIDATE{
        long long resultvecX[8];
        _mm512_storeu_si512((void *) resultvecX, resultvec);
        //    printf("%ld %ld %ld %ld %ld %ld %ld %ld\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
        //    printf("%lx %lx %lx %lx %lx %lx %lx %lx\n", resultvecX[0],resultvecX[1],resultvecX[2],resultvecX[3],resultvecX[4],resultvecX[5],resultvecX[6],resultvecX[7]);
        result += resultvecX[0] + resultvecX[1] + resultvecX[2] + resultvecX[3] + resultvecX[4] + resultvecX[5] + resultvecX[6] + resultvecX[7];
//        printf("par = %d: %d %d\n", par, result, _mm512_reduce_add_epi64(resultvec));
      }
    EGO_END
    printf("result %d\n", result);
}

int test_ego_loop2(int n, int d, double epsilon, double *array){
    long long result = 0;
    EGO_PARALLEL_TRAN(n, d, epsilon, 5, array)
        #pragma omp parallel for reduction(+:result)
    EGO_PREPARE
        veci64 resultvec = _mm512_setzero_si512();
    EGO_LOOP_TRAN{
        resultvec += _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63))
                +  _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63))
                +  _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63))
                +  _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63))
                +  _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63))
                +  _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63))
                +  _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63))
                +  _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63)) ;
    }
    EGO_CONSOLIDATE{
        result += _mm512_reduce_add_epi64(resultvec);
//        double testres[8] __attribute__((aligned(64)));
//        _mm512_store_epi64(testres, resultvec);
//        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
    }
    EGO_END
    printf("result %d\n", result);
}

int test_ego_loop3(int n, int d, double epsilon, double *array){
    size_t result = 0;
    EGO_PARALLEL_TRAN(n, d, epsilon, 14, array)
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
    printf("result %zu\n", result);
}

//void test(int i, int j, double * array){
//    vec sum0=_mm512_setzero_pd ;
//    vec sum1=_mm512_setzero_pd ;
//    vec sum2=_mm512_setzero_pd ;
//    vec sum3=_mm512_setzero_pd ;
//    vec sum4=_mm512_setzero_pd ;
//    vec sum5=_mm512_setzero_pd ;
//    vec sum6=_mm512_setzero_pd ;
//    vec sum7=_mm512_setzero_pd ;
//    veci64 increment = _mm512_set1_epi64(1ll) ;
//    for(int k=0 ; k<d ; k++){
//        vec vi=_mm512_load_pd(array+k*n+i*8);
//        vec vj=_mm512_load_pd(array+k*n+j*8);
//        sum0 += vi * _mm512_permutexvar_pd(increment, vj);
//        veci64 cur = increment + increment ;
//        sum1 += vi * _mm512_permutexvar_pd(cur, vj);
//        cur += increment ;
//        sum2 += vi * _mm512_permutexvar_pd(cur, vj);
//        cur += increment ;
//        sum3 += vi * _mm512_permutexvar_pd(cur, vj);
//        cur += increment ;
//        sum4 += vi * _mm512_permutexvar_pd(cur, vj);
//        cur += increment ;
//        sum5 += vi * _mm512_permutexvar_pd(cur, vj);
//        cur += increment ;
//        sum6 += vi * _mm512_permutexvar_pd(cur, vj);
//        cur += increment ;
//        sum7 += vi * _mm512_permutexvar_pd(cur, vj);
//
//    }
//    vec via = _mm512_load_pd(self + i * 40);
//    vec vib = _mm512_load_pd(self + i * 40 + 8);
//    vec vic = _mm512_load_pd(self + i * 40 + 16);
//    vec vid = _mm512_load_pd(self + i * 40 + 24);
//    vec vie = _mm512_load_pd(self + i * 40 + 32);
//    vec vj = _mm512_set1_pd(self + j * 40);
//    vec sum0a = vj+via; vec sum0b = vj+vib ; vec sum0c = vj+vic ; vec sum0d = vj+vid ; vec sum0e = vj+vie ;
//    vec vj = _mm512_set1_pd(self + j * 40 + 1);
//    vec sum1a = vj+via; vec sum1b = vj+vib ; vec sum1c = vj+vic ; vec sum1d = vj+vid ; vec sum1e = vj+vie ;
//    vec vj = _mm512_set1_pd(self + j * 40 + 2);
//    vec sum2a = vj+via; vec sum2b = vj+vib ; vec sum2c = vj+vic ; vec sum2d = vj+vid ; vec sum2e = vj+vie ;
//    vec vj = _mm512_set1_pd(self + j * 40 + 3);
//    vec sum3a = vj+via; vec sum3b = vj+vib ; vec sum3c = vj+vic ; vec sum3d = vj+vid ; vec sum3e = vj+vie ;
//    vec vj = _mm512_set1_pd(self + j * 40 + 4);
//    vec sum4a = vj+via; vec sum4b = vj+vib ; vec sum4c = vj+vic ; vec sum4d = vj+vid ; vec sum4e = vj+vie ;
//    for (int k=0 ; k<d ; k++){
//        via = _mm512_load_pd(array+(i*d+k)*40);
//        vib = _mm512_load_pd(array+(i*d+k)*40+8);
//        vic = _mm512_load_pd(array+(i*d+k)*40+16);
//        vid = _mm512_load_pd(array+(i*d+k)*40+24);
//        vie = _mm512_load_pd(array+(i*d+k)*40+32);
//        vj = _mm512_set1_pd(array+(j*d+k)*40);
//        sum0a += vj*via ; sum0b += vj*vib ; sum0c += vj*vic ; sum0d += vj*vid ; sum0e += vj*vid ;
//        vj = _mm512_set1_pd(array+(j*d+k)*40+1);
//        sum1a += vj*via ; sum1b += vj*vib ; sum1c += vj*vic ; sum1d += vj*vid ; sum1e += vj*vid ;
//        vj = _mm512_set1_pd(array+(j*d+k)*40+2);
//        sum2a += vj*via ; sum2b += vj*vib ; sum2c += vj*vic ; sum2d += vj*vid ; sum2e += vj*vid ;
//        vj = _mm512_set1_pd(array+(j*d+k)*40+3);
//        sum3a += vj*via ; sum3b += vj*vib ; sum3c += vj*vic ; sum3d += vj*vid ; sum3e += vj*vid ;
//        vj = _mm512_set1_pd(array+(j*d+k)*40+4);
//        sum4a += vj*via ; sum4b += vj*vib ; sum4c += vj*vic ; sum4d += vj*vid ; sum4e += vj*vid ;
//
//}

void transClosureDens(int n, int d, double epsilon, int MinPts, double *array, long long *counter){
    listArray = (ListEl **) malloc(n * sizeof (ListEl*));
    listArray[0] = (ListEl *) malloc(n * sizeof(ListEl));
    for (int i=0 ; i<n ; i++){
        listArray[i] = listArray[0] + i;
        listArray[i]->n=1;
        listArray[i]->minid=i ;
        listArray[i]->left = (ListEl*)0;
        omp_init_lock(&listArray[i]->lock);
    }
    EGO_PARALLEL(n, d, epsilon, array)
        #pragma omp parallel for
    EGO_PREPARE
        double * temp = mallocA64(8*sizeof(double));
    EGO_LOOP{
        _mm512_store_pd(temp, indicator);
        for(int k=0; k<8 ; k++)
            if(temp[k] >= 0) {
                if(counter[i*8+k] >= MinPts){
                    if(counter[j*8+k] >= MinPts)
                        unifyList(i*8+k, j*8+k, 0);
                    else if(listArray[j*8+k]->n == 1) // this must be protected by a lock
                        unifyList(i*8+k, j*8+k, 0);
                } else {
                    if(counter[j*8+k] >= MinPts)
                        if(listArray[j*8+k]->n == 1) // this must be protected by a lock
                            unifyList(i*8+k, j*8+k, 0);
                }
            }
    }
    EGO_CONSOLIDATE{
    }
    EGO_END
}

long long * corePoints(int n, int d, double epsilon, double *array){
    long long * counter = (long long*)callocA64(n*sizeof(long long));
    omp_lock_t locks[2*NUM_THREADS] ;
    for (int i=0 ; i<2*NUM_THREADS ; i++)
        omp_init_lock(locks+i);
    EGO_PARALLEL(n, d, epsilon, array)
        #pragma omp parallel for
    EGO_PREPARE
        long long * locct = (long long*)callocA64(n*sizeof(long long));
        long long * temp = (long long*)mallocA64(64*sizeof(long long));
    EGO_LOOP{
        //if(i==0)print_vector(indicator);
        veci64 vi8 = _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(indicator), 63));
        _mm512_store_si512(locct + j*8, _mm512_add_epi64(_mm512_load_si512(locct + j*8), vi8)) ;
        if(((III+II)&7) == 7){
            veci64 vi1 = _mm512_load_si512(temp) ;
            veci64 vi2 = _mm512_load_si512(temp+8) ;
            veci64 vi3 = _mm512_load_si512(temp+16) ;
            veci64 vi4 = _mm512_load_si512(temp+24) ;
            veci64 vi5 = _mm512_load_si512(temp+32) ;
            veci64 vi6 = _mm512_load_si512(temp+40) ;
            veci64 vi7 = _mm512_load_si512(temp+48) ;
            transposeAVX512i(vi1,vi2,vi3,vi4,vi5,vi6,vi7,vi8);
            _mm512_store_si512(locct + i*8, _mm512_add_epi64(_mm512_load_si512(locct + i*8), vi1+vi2+vi3+vi4+vi5+vi6+vi7+vi8)) ;
        } else {
            _mm512_store_si512(temp + ((III+II)&7)*8, vi8);
        }
    }
    EGO_CONSOLIDATE{
        int segments_added=0;
        int curseg = 2*(omp_get_thread_num()%NUM_THREADS) ;
        char * segment=(char*) calloc(2*NUM_THREADS,1);
        while(segments_added < 2 * NUM_THREADS){
            if(!segment[curseg] && omp_test_lock(locks+curseg)){
                for(int i=curseg*n/NUM_THREADS/16*8 ; i<(curseg+1)*n/NUM_THREADS/16*8 ; i+=8)
                    _mm512_store_si512(counter+i, _mm512_add_epi64(_mm512_load_si512(counter+i),_mm512_load_si512(locct+i)));
                omp_unset_lock(locks+curseg);
                segments_added++;
                segment[curseg] = 1;
            }
            if(++curseg >= 2*NUM_THREADS)
                curseg = 0;
        }
    }
    EGO_END
    return counter;
}
long long * corePoints2(int n, int d, double epsilon, double *array){
    long long * counter = (long long*)callocA64(n*sizeof(long long));
    omp_lock_t locks[2*NUM_THREADS] ;
    for (int i=0 ; i<2*NUM_THREADS ; i++)
        omp_init_lock(locks+i);
    EGO_PARALLEL_TRAN(n, d, epsilon, 2, array)
        #pragma omp parallel for
    EGO_PREPARE
        long long * locct = (long long*)callocA64(n*sizeof(long long));
        long long * temp = (long long*)mallocA64(64*sizeof(long long));
        veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};
        veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};
        veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};
        veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};
        veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};
        veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};
    EGO_LOOP_TRAN{
        veci64 vi1 = _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63));
        veci64 vi2 = _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63));
        veci64 vi3 = _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63));
        veci64 vi4 = _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63));
        veci64 vi5 = _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63));
        veci64 vi6 = _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63));
        veci64 vi7 = _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63));
        veci64 vi8 = _mm512_xor_epi64(_mm512_set1_epi64(1ll), _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63));
        _mm512_store_si512(locct + j*8, _mm512_add_epi64(_mm512_load_si512(locct + j*8), vi1+vi2+vi3+vi4+vi5+vi6+vi7+vi8)) ;
        transposeAVX512i(vi1,vi2,vi3,vi4,vi5,vi6,vi7,vi8);
        _mm512_store_si512(locct + i*8, _mm512_add_epi64(_mm512_load_si512(locct + i*8), vi1+vi2+vi3+vi4+vi5+vi6+vi7+vi8)) ;
    }
    EGO_CONSOLIDATE{
        int segments_added=0;
        int curseg = 2*(omp_get_thread_num()%NUM_THREADS) ;
        char * segment=(char*) calloc(2*NUM_THREADS,1);
        while(segments_added < 2 * NUM_THREADS){
            if(!segment[curseg] && omp_test_lock(locks+curseg)){
                for(int i=curseg*n/NUM_THREADS/16*8 ; i<(curseg+1)*n/NUM_THREADS/16*8 ; i+=8)
                    _mm512_store_si512(counter+i, _mm512_add_epi64(_mm512_load_si512(counter+i),_mm512_load_si512(locct+i)));
                omp_unset_lock(locks+curseg);
                segments_added++;
                segment[curseg] = 1;
            }
            if(++curseg >= 2*NUM_THREADS)
                curseg = 0;
        }
    }
    EGO_END
    return counter;
}
long long * corePoints3(int n, int d, double epsilon, double *array){
    long long * counter = (long long*)callocA64(n*sizeof(long long));
    omp_lock_t locks[2*NUM_THREADS] ;
    for (int i=0 ; i<2*NUM_THREADS ; i++)
        omp_init_lock(locks+i);
    // Aufruf EGO-Join EGO_PARALLEL_TRAN ... EGO_END_TRAN
    EGO_PARALLEL_TRAN(n, d, epsilon, 2, array)
        #pragma omp parallel for
    EGO_PREPARE
        long long * locct = (long long*)callocA64(n*sizeof(long long));
        long long * temp = (long long*)mallocA64(64*sizeof(long long));
        veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};
        veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};
        veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};
        veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};
        veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};
        veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};
        veci64 eights = _mm512_set1_epi64(8ll) ;
    EGO_LOOP_TRAN{
        veci64 vi1 = _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63);
        veci64 vi2 = _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63);
        veci64 vi3 = _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63);
        veci64 vi4 = _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63);
        veci64 vi5 = _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63);
        veci64 vi6 = _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63);
        veci64 vi7 = _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63);
        veci64 vi8 = _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63);
        _mm512_store_si512(locct + i*8, _mm512_load_si512(locct + i*8) + eights - (vi1+vi2+vi3+vi4+vi5+vi6+vi7+vi8)) ;
        transposeAVX512i(vi1,vi2,vi3,vi4,vi5,vi6,vi7,vi8);
        _mm512_store_si512(locct + j*8, _mm512_load_si512(locct + j*8) + eights - (vi1+vi2+vi3+vi4+vi5+vi6+vi7+vi8)) ;
    }
    EGO_CONSOLIDATE{
        int segments_added=0;
        int curseg = 2*(omp_get_thread_num()%NUM_THREADS) ;
        char * segment=(char*) calloc(2*NUM_THREADS,1);
        while(segments_added < 2 * NUM_THREADS){
            if(!segment[curseg] && omp_test_lock(locks+curseg)){
                for(int i=curseg*n/NUM_THREADS/16*8 ; i<(curseg+1)*n/NUM_THREADS/16*8 ; i+=8)
                    _mm512_store_si512(counter+i, _mm512_add_epi64(_mm512_load_si512(counter+i),_mm512_load_si512(locct+i)));
                omp_unset_lock(locks+curseg);
                segments_added++;
                segment[curseg] = 1;
            }
            if(++curseg >= 2*NUM_THREADS)
                curseg = 0;
        }
    }
    EGO_END_TRAN
    return counter;
}

inline void bla(int i, int j, int MinPts, long long *counter, vec indicator) {
    double temp[8]__attribute__((aligned(64)));
    _mm512_store_pd(temp, indicator);
    for (int k = 0; k < 8; k++)
        if (temp[k] >= 0) {
            unifyListDbscan(i+k, j, counter[i+k] >= MinPts, counter[j] >= MinPts);
        }
}
void dbscan(int n, int d, double epsilon, int MinPts, double *array){
    long long * counter = (long long*)callocA64(n*sizeof(long long));
//    unifyCounter = 0;
//    histo = (int *) calloc(HISTO_SIZE, sizeof(int));
    omp_lock_t locks[2*NUM_THREADS] ;
    for (int i=0 ; i<2*NUM_THREADS ; i++)
        omp_init_lock(locks+i);
    EGO_PARALLEL_TRAN(n, d, epsilon, 5, array)
        #pragma omp parallel for
    EGO_PREPARE
        long long * locct = (long long*)callocA64(n*sizeof(long long));
        //long long * temp = (long long*)mallocA64(64*sizeof(long long));
        veci64 mask1 = {0, 8, 2, 10, 4, 12, 6, 14};
        veci64 mask2 = {9, 1, 11, 3, 13, 5, 15, 7};
        veci64 mask3 = {0, 1, 8, 9, 4, 5, 12, 13};
        veci64 mask4 = {10, 11, 2, 3, 14, 15, 6, 7};
        veci64 mask5 = {0, 1, 2, 3, 8, 9, 10, 11};
        veci64 mask6 = {12, 13, 14, 15, 4, 5, 6, 7};
        veci64 eights = _mm512_set1_epi64(8ll) ;
    EGO_LOOP_TRAN{
        veci64 vi1 = _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63);
        veci64 vi2 = _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63);
        veci64 vi3 = _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63);
        veci64 vi4 = _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63);
        veci64 vi5 = _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63);
        veci64 vi6 = _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63);
        veci64 vi7 = _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63);
        veci64 vi8 = _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63);
        _mm512_store_si512(locct + i*8, _mm512_load_si512(locct + i*8) + eights - (vi1+vi2+vi3+vi4+vi5+vi6+vi7+vi8)) ;
        transposeAVX512i(vi1,vi2,vi3,vi4,vi5,vi6,vi7,vi8);
        _mm512_store_si512(locct + j*8, _mm512_load_si512(locct + j*8) + eights - (vi1+vi2+vi3+vi4+vi5+vi6+vi7+vi8)) ;
    }
    EGO_CONSOLIDATE{
        int segments_added=0;
        int curseg = 2*(omp_get_thread_num()%NUM_THREADS) ;
        char * segment=(char*) calloc(2*NUM_THREADS,1);
        while(segments_added < 2 * NUM_THREADS){
            if(!segment[curseg] && omp_test_lock(locks+curseg)){
                for(int i=curseg*n/NUM_THREADS/16*8 ; i<(curseg+1)*n/NUM_THREADS/16*8 ; i+=8)
                    _mm512_store_si512(counter+i, _mm512_add_epi64(_mm512_load_si512(counter+i),_mm512_load_si512(locct+i)));
                omp_unset_lock(locks+curseg);
                segments_added++;
                segment[curseg] = 1;
            }
            if(++curseg >= 2*NUM_THREADS)
                curseg = 0;
        }
    }
    EGO_END_RESUME
        printf("End count points\n");
        for (int i=0 ; i<20 ; i++)printf("%ld ", counter[i]); printf("\n");
        for (int i=n-20; i<n ; i++)printf("%ld ", counter[i]); printf("\n");
        long long allall=0; for(int i=0 ; i<n ; i++)allall+=counter[i]; printf("allall=%ld avg=%ld\n", allall, allall/n);
        listArray = (ListEl **) malloc(n * sizeof (ListEl*));
        listArray[0] = (ListEl *) malloc(n * sizeof(ListEl));
        for (int i=0 ; i<n ; i++){
            listArray[i] = listArray[0] + i;
            listArray[i]->n=1;
            listArray[i]->minid=i ;
            listArray[i]->left = (ListEl*)0;
            omp_init_lock(&listArray[i]->lock);
        }
    EGO_PARALLEL_TRAN_RESUME
        #pragma omp parallel for
    EGO_PREPARE
    EGO_LOOP_TRAN{
//        if(i%256 == 0 && j%256 == 0) printf("%d %d %d\n", omp_get_thread_num(),i,j);
        bla(i*8,j*8,MinPts,counter,sum1);
        bla(i*8,j*8+1,MinPts,counter,sum2);
        bla(i*8,j*8+2,MinPts,counter,sum3);
        bla(i*8,j*8+3,MinPts,counter,sum4);
        bla(i*8,j*8+4,MinPts,counter,sum5);
        bla(i*8,j*8+5,MinPts,counter,sum6);
        bla(i*8,j*8+6,MinPts,counter,sum7);
        bla(i*8,j*8+7,MinPts,counter,sum8);
    }
    EGO_CONSOLIDATE{
    }
    EGO_END_TRAN
//    printf("UnifyCounter = %ld\n", unifyCounter) ;
//    for(int i=0 ; i<HISTO_SIZE ; i++)
//        printf("%d ", histo[i]);
//    printf("\n");
//    for(int i=0 ; i<20 ; i++)
//        printf("%d %d %d %d\n", i, listArray[i]->minid, listArray[i]->n, counter[i]) ;
//    printf("\n");
    for(int i=0 ; i<n ; i++)
        if(listArray[i]->minid == i && listArray[i]->n > 1000){
            int ct = 0;
            for (int j=i ; j<n ; j++)
                if(listArray[j]->minid == i && counter[j] >= MinPts) ct++ ;
            printf("%d %d(core %d) %d ", listArray[i]->minid, listArray[i]->n, ct, counter[i]) ;
            printf("(%d %d %d ...)\n", listArray[i]->left->minid, listArray[i]->left->left->minid, listArray[i]->left->left->left->minid);
        }
//    int ct=0 ;
//    for (int i=0 ; i<n && ct < 200 ; i++) if(listArray[i]->minid == 19996){
//        ct ++;
//        for(int j=0 ; j<d ; j++)
//            printf("%9.6lf ", array[(i/8*d+j)*8+i%8]);
//        printf("\n");
//    }
//    printf("\n");
//    ct=0 ;
//    for (int i=0 ; i<n && ct < 200 ; i++) if(listArray[i]->minid == 0){
//        ct ++;
//        for(int j=0 ; j<d ; j++)
//            printf("%9.6lf ", array[(i/8*d+j)*8+i%8]);
//        printf("\n");
//    }
}
#ifndef XFOLD
#define XFOLD 1
#endif

void testInt(int m){
    int neun = 3;
    while (neun < m)
        neun*=3 ;
    neun /= 3 ;
    for (int i=0 ; i<neun; i++){
        int mi = ((i+1)*2*m+neun)/neun/2 - (i*2*m+neun)/neun/2;
        for (int j=0 ; j<mi ; j++)
            printf("%d ", 3*i+(j*2*3+mi)/mi/2);
        printf("\n");
    }
    printf("\n");
}
void testInt2(int m, int verbose){
    int neun = 1;
    while (neun < m)
        neun*=3 ;
    int lasti = -1;
    int testarray[neun];
    for (int k=0 ; k<m ; k++){
        if (verbose && k==m/2)
            printf("----------------------\n");
        int i=(k*neun*2+m)/m/2;
        int i2=((k-1)*neun*2+m)/m/2;
        if((i+3)/3 > (i2+3)/3)
            i-= i%3 ;
        if(i==lasti)
            printf("ERROR! m=%d %d %d %d %d\n", m, k, (k*neun*2+m)/m/2, i2, i);
        lasti=i;
        testarray[i] = m;
        int umrechnung = 0;
        int power = 1;
        for (int j=i; j>0 ; j/=3){
            umrechnung += power*(j%3) ;
            power *=10;
        }
        if(verbose)printf(m>27?"%3d %3d %3d %3d   %04d\n":m>9?"%3d %3d %3d %3d   %03d\n":"%3d %3d %3d %3d   %02d\n", k, (k*neun*2+m)/m/2, i2, i, umrechnung);
    }
    if(verbose)printf("\n");
    for(int i=0 ; i<neun ; i+=3)
        if(testarray[i]!=m)
            printf("ERROR2! m=%d i=%d\n", m, i);
}
// void generateSinData(int n, int d, int k, double noise, double* array){
//     double rot[k][2][d];
//     double trans[k][d];
//     for (int i=0 ; i<k*2*d ; i++)
//         rot[0][0][i] = drand48() + drand48() + drand48() - 1.5;
//     for(int i=0 ; i<k ; i++){
//         double h=0;
//         for (int j=0 ; j<d ; j++) h+=rot[i][0][j]*rot[i][0][j];
//         h = sqrt(h) ;
//         for (int j=0 ; j<d ; j++) rot[i][0][j]/=h;
//         h=0;
//         for (int j=0 ; j<d ; j++) h += rot[i][0][j]*rot[i][1][j];
//         for (int j=0 ; j<d ; j++) rot[i][1][j] -= h * rot[i][0][j];
//         h=0;
//         for (int j=0 ; j<d ; j++) h+=rot[i][1][j]*rot[i][1][j];
//         h = sqrt(h) ;
//         for (int j=0 ; j<d ; j++) rot[i][1][j]/=h;
//     }
//     for (int i=0 ; i<k*d ; i++)
//         trans[0][i] = (drand48() + drand48() + drand48() - 1.5) * 3.;
//     for(int i=0 ; i<n ; i++){
//         if(drand48()<noise){
//             for (int j=0 ; j<d ; j++)
//                 array[i*d+j] = drand48() + drand48() + drand48() - 1.5;
//         } else {
//             int c=drand48()*k;
//             double x=drand48()*2.0-1.0;
//             double y=sin(3.14*x);
//             for (int j=0 ; j<d ; j++)
//                 array[i*d+j] = rot[c][0][j] * x + rot[c][1][j] * y + trans[c][j] + drand48() / 3 / d;
//         }
// //        for(int j=0 ; j<d ; j++)
// //            printf("%6.3f ", array[i*d+j]);
// //        printf("\n");
//     }
// }
int * dbscan_processed;
int dbscan_n;
int dbscan_d;
//double dbscan_epsilon;
double dbscan_eps_sq;
int dbscan_minPts;
double * dbscan_array;
int dbscan_curId;
int dbscan_cores;
int dbscan_stripes;
int **dbscan_lower;
int **dbscan_upper;

inline int atomic_mark_processed(int p){
    int result;
#pragma omp critical
    {
        result = dbscan_processed[p];
        if(!result)
            dbscan_processed[p] = dbscan_curId;
    }
    return result;
}
void expandCluster(int p) {
#pragma omp task
    {
//        if (!dbscan_processed[p]) {
//            dbscan_processed[p] = dbscan_curId;
        if(!atomic_mark_processed(p)){
            int counter = 0;
            int * memory = (int*)malloc(dbscan_minPts * sizeof(int));
//            int memory[dbscan_minPts];
//            for (int i = 0; i < dbscan_n; i++)
            for (int s = 0; s < dbscan_stripes; s++)
                for (int i = dbscan_lower[s][p]; i < dbscan_upper[s][p]; i++)
                    if (!dbscan_processed[i] || counter < dbscan_minPts) {
                        double distance = 0;
                        for (int j = 0; j < dbscan_d && distance <= dbscan_eps_sq ; j++) {
                            double h = dbscan_array[i * dbscan_d + j] - dbscan_array[p * dbscan_d + j];
                            distance += h*h;
                        }
                        if (distance <= dbscan_eps_sq) {
                            {
                                if (counter < dbscan_minPts) {
                                    memory[counter] = i;
                                } else {
                                    expandCluster(i);
                                }
                                counter++;
                            }
                        }
                    }
            if (counter >= dbscan_minPts) {
#pragma omp critical(countcores)
                {
                dbscan_cores++;
                }
                for (int j = 0; j < dbscan_minPts; j++)
                    if (!dbscan_processed[memory[j]]) {
                        expandCluster(memory[j]);
                    }
            }
            free(memory);
        }
    }
}

void dbscan_simple(int n, int d, double epsilon, int minPts, double * array) {
    dbscan_stripes = 9;
    dbscan_n = n;
    dbscan_d = d;
    dbscan_eps_sq = epsilon*epsilon;
    dbscan_minPts = minPts;
    dbscan_array = array;
    dbscan_curId = 1;
    dbscan_processed = (int*) calloc(n * sizeof (int), 1);
    epsilonGridOrdering(n, d, epsilon, array);
    dbscan_lower = (int **) malloc(dbscan_stripes * sizeof (int*));
    dbscan_upper = (int **) malloc(dbscan_stripes * sizeof (int*));
    prepareStripesFull(n, d, dbscan_stripes, epsilon, array, dbscan_lower, dbscan_upper);
    int localcounter[NUM_THREADS];
    int localalloc[NUM_THREADS];
    int *localmem[NUM_THREADS];
    for (int par = 0; par < NUM_THREADS; par++) {
        localalloc[par] = 100;
        localmem[par] = (int *) malloc(localalloc[par] * sizeof (int));
    }
    for (int p = 0; p < n; p++) {
        dbscan_cores = 0;
        if (!dbscan_processed[p]) {
            dbscan_processed[p] = dbscan_curId;
            int allload = 0;
            for (int s = 0; s < dbscan_stripes; s++)
                allload += dbscan_upper[s][p] - dbscan_lower[s][p];
            int counter = 0;
//            if(p%65536 == 0)
//                printf("%d allload=%d\n", p, allload);
#pragma omp parallel for reduction(+:counter)
            for (int par = 0; par < NUM_THREADS; par++) {
                int mycumload = 0;
                localcounter[par] = 0;
                for (int s = 0; s < dbscan_stripes; s++) {
                    if (mycumload >= (par + 1) * allload / NUM_THREADS)
                        break;
                    if (mycumload + dbscan_upper[s][p] - dbscan_lower[s][p] < par * allload / NUM_THREADS) {
                        mycumload += dbscan_upper[s][p] - dbscan_lower[s][p];
                        continue;
                    }
                    int firstadd = par * allload / NUM_THREADS - mycumload;
                    if (firstadd < 0) firstadd = 0;
                    mycumload += firstadd;
                    for (int i = dbscan_lower[s][p] + firstadd;
                            i < dbscan_upper[s][p] && mycumload < (par + 1) * allload / NUM_THREADS;
                            i++, mycumload++) {
                        double distance = 0;
                        for (int j = 0; j < dbscan_d && distance <= dbscan_eps_sq; j++) {
                            double h = dbscan_array[i * dbscan_d + j] - dbscan_array[p * dbscan_d + j];
                            distance += h*h;
                        }
                        if (distance <= dbscan_eps_sq) {
                            localcounter[par]++;
                            if (localcounter[par] > localalloc[par]) {
                                localalloc[par] *= 2;
                                localmem[par] = (int *) realloc(localmem[par], localalloc[par] * sizeof (int));
                            }
                            localmem[par][localcounter[par] - 1] = i;
                        }
                    }
                }
                counter += localcounter[par];
            }
            if (counter < dbscan_minPts) {
                dbscan_processed[p] = 0; //will be re-visited at most once!
            } else {
#pragma omp critical(countcores)
                {
                    dbscan_cores++;
                }
//                printf("Start Cluster %d\n", p);
#pragma omp parallel
#pragma omp master
                for (int par = 0; par < NUM_THREADS; par++) {
                    for (int i = 0; i < localcounter[par]; i++)
                        if (!dbscan_processed[localmem[par][i]]) {
                            expandCluster(localmem[par][i]);
                        }
                }
#pragma omp taskwait
            }
        }
        if (dbscan_cores) {
            if (dbscan_cores > 500)
                printf("Cluster %d p=%d crs=%d\n", dbscan_curId, p, dbscan_cores);
            dbscan_curId++;
        }
    }
    for (int par = 0; par < NUM_THREADS; par++){
//        printf("%d ",localalloc[par]);
        free(localmem[par]);
    }
}

void dbscan_semaphor(int n, int d, double epsilon, int minPts, double * array) {
    int stripes = 9;
    epsilonGridOrdering(n, d, epsilon, array);
    int **lower = (int **) malloc(stripes * sizeof (int*));
    int **upper = (int **) malloc(stripes * sizeof (int*));
    prepareStripesFull(n, d, stripes, epsilon, array, lower, upper);
    int * processed = (int *) calloc(n, sizeof (int));
    int * seedlist = (int *) malloc(n * sizeof (int));
    double eps_sq = epsilon * epsilon ;
    int curId = 1;
    for (int p = 0; p < n; p++) {
        int cores = 0;
        if (!processed[p]) {
            processed[p] = curId;
            int allload = 0;
            for (int s = 0; s < stripes; s++)
                allload += upper[s][p] - lower[s][p];
            int seedcounter = 0;
#pragma omp parallel for
            for (int par = 0; par < NUM_THREADS; par++) {
                int mycumload = 0;
                for (int s = 0; s < stripes; s++) {
                    if (mycumload >= (par + 1) * allload / NUM_THREADS)
                        break;
                    if (mycumload + upper[s][p] - lower[s][p] < par * allload / NUM_THREADS) {
                        mycumload += upper[s][p] - lower[s][p];
                        continue;
                    }
                    int firstadd = par * allload / NUM_THREADS - mycumload;
                    if (firstadd < 0) firstadd = 0;
                    mycumload += firstadd;
                    for (int i = lower[s][p] + firstadd;
                            i < upper[s][p] && mycumload < (par + 1) * allload / NUM_THREADS;
                            i++, mycumload++) {
                        double distance = 0;
                        for (int j = 0; j < d && distance <= eps_sq; j++) {
                            double h = array[i * d + j] - array[p * d + j];
                            distance += h*h;
                        }
                        if (distance <= eps_sq) {
#pragma omp critical
                            {
                                seedlist[seedcounter++] = i;
                            }
                        }
                    }
                }
            }
            if (seedcounter >= minPts) {
                processed[p] = curId;
                omp_lock_t delay;
                omp_init_lock(&delay);
                int running = 0;
                int die = 0;
                cores = 1;
                printf("Start Cluster %d (%d init %d obj)\n", curId, p, seedcounter);
#pragma omp parallel for
                for (int par = 0; par < NUM_THREADS; par++) {
                    int ii;
                    while (!die) {
                        // WAIT for the counting semaphor
                        omp_set_lock(&delay);
#pragma omp critical(mutex)
                        {
                            if (seedcounter + running <= 0) {
                                die = 1;
                            }
                            if (!die) {
                                ii = seedlist[--seedcounter];
                                running++;
                                if (seedcounter > 0)
                                    omp_unset_lock(&delay);
                            }
//                            printf("PID %d: seedcounter: %d running=%d obj %d\n", omp_get_thread_num(), seedcounter, running, ii);
                        }
                        if (die) break;
                        int memory[minPts];
                        int counter = 0;
                        //            for (int i = 0; i < n; i++)
                        for (int s = 0; s < stripes; s++)
                            for (int i = lower[s][ii]; i < upper[s][ii]; i++)
                                if (!processed[i] || counter < minPts) {
                                    double distance = 0;
                                    for (int j = 0; j < d && distance <= eps_sq; j++) {
                                        double h = array[i * d + j] - array[ii * d + j];
                                        distance += h*h;
                                    }
                                    if (distance <= eps_sq) {
                                        if (counter < minPts) {
                                            memory[counter] = i;
                                        } else {
                                            // SIGNAL the counting semaphor
#pragma omp critical(mutex)
                                            {
                                                if (!processed[i]) {
                                                    processed[i] = curId;
                                                    seedlist[seedcounter++] = i;
                                                    if (seedcounter == 1) omp_unset_lock(&delay);
                                                }
                                            }
                                        }
                                        counter++;
                                    }
                                }
#pragma omp critical(mutex)
                        {
                            if (counter >= minPts) {
                                // SIGNAL the counting semaphor
                                cores++;
                                for (int i = 0; i < minPts; i++)
                                    if (!processed[memory[i]]) {
                                        processed[memory[i]] = curId;
                                        seedlist[seedcounter++] = memory[i];
                                        if (seedcounter == 1) omp_unset_lock(&delay);
                                    }
                            }
                            running--;
                            if (seedcounter + running <= 0) {
                                die = 1;
                            }
//                            printf("PID %d: seedcounter: %d running=%d obj %d END\n", omp_get_thread_num(), seedcounter, running, ii);
                        }
                    }
                    omp_unset_lock(&delay);
//                    printf("Died %d running=%d\n", omp_get_thread_num(), running);
                }
                curId++;
                if (cores > 500)
                    printf("Cluster %d p=%d crs=%d\n", curId, p, cores);
            }
        }
    }
}
int main(int argc, char** argv) {
//    generateSinData(300,3,3,0.05,(double*)malloc(300*3*sizeof(double)));
//    exit(0);
    int n = 11620300; // durch 8 teilbar für test_ego_loop
    int d = 57;
    double epsilon = 0.06;
    omp_set_num_threads(NUM_THREADS);
    printf("n=%d d=%d eps=%lf threads=%d omp_threads=%d\n",n,d,epsilon,NUM_THREADS,omp_get_max_threads());
    double * array = (double*) mallocA64((n+7)/8*8 * sizeof (double) * d + 16384);
//    generateSinData(n,d,5,0.001,array); FILE *ff = fopen("sinclust4-11-5", "w"); fwrite(array, n*d, sizeof(double), ff); fclose(ff);
//    exit(0);
    FILE * f3 = fopen("/home/share/bigcross/bigCross_11620300x57_normalized_reorder_0.06.bin", "r");
//    FILE * f3 = fopen("sinclust283", "r");
//    FILE * f3 = fopen("bmatrix", "r");
    fread(array, n*d, sizeof (double), f3);
    fclose(f3);
//    for(int i=n-16 ; i<n+2 ; i++){
//        for (int j=0 ; j<d ; j++)
//            printf("%8.5f ", array[i*d+j]);
//        printf("\n");
//    }
//    stop(); stopc();
//    transpose_8xd(n,d,array);
//    printf("transpose_8xd| %8.3f | %8.3f\n", stopc(), stop());
//    for(int i=0 ; i<2*d ; i++){
//        for (int j=0 ; j<8 ; j++)
//            printf("%8.5f ", array[i*8+j]);
//        printf("\n");
//    }
//    stop(); stopc();
//    transpose_dx8(n,d,array);
//    printf("transpose_dx8| %8.3f | %8.3f\n", stopc(), stop());
//    for(int i=n-16 ; i<n+2 ; i++){
//        for (int j=0 ; j<d ; j++)
//            printf("%8.5f ", array[i*d+j]);
//        printf("\n");
//    }
//    exit(0);

//    stop(); stopc();
//    corePoints3(n,d,epsilon,array);
//    printf("corePoints3 | %8.3f | %8.3f\n", stopc(), stop());

//     printf("vorher Test Hilbert Loop\n");
//     stop();
//     stopc();
//        int i=0 ; int j=0;
// //       FUR_HILBERT_START(i,j,0,300000,0,300000){
//        FGF_HILBERT_FOR(i, j, 1 << 30 , 1 << 30, i >=1000000000 && i < 1000000010 && j >= 1000000000 && j < 1000000010, FURHIL_lb0 <=1000000010 && FURHIL_ub0 >=1000000000 && FURHIL_lb1<=1000000010 && FURHIL_ub1 >= 1000000000) {
// //           if(i >=200000 && i < 200010 && j >= 200000 && j < 200010 || HILLOOP_hilbert % 10000000ULL == 0ULL)
// //           if(HILLOOP_hilbert == 200000000ULL)
// //           if(HILLOOP_hilbert % 10000000ULL == 0ULL)
//            printf("%ld %d %d \n", HILLOOP_hilbert, i,j);
// //           printf("%ld %d %d %d %d %d %d %d %d\n", HILLOOP_hilbert, i,j, HILLOOP_isize, HILLOOP_jsize, (HILLOOP_J + 1) * (HILLOOP_jcur-HILLOOP_j57), HILLOOP_j57, HILLOOP_J, HILLOOP_jcur);
// //           if(HILLOOP_hilbert == 0ULL){
// //               printf("%ld\n", HILLOOP_stop);
//
//
//        }FGF_HILBERT_END(i,j)
// //       }FUR_HILBERT_END(i,j)
//     printf("FURHIL | %8.3f | %8.3f\n", stopc(), stop());

    stop();
    stopc();
    // dbscan(n,d,epsilon,10,array);
    // printf("dbscan| %8.3f | %8.3f\n", stopc(), stop());
    // dbscan_semaphor(n,d,epsilon,11,array);
    // printf("dbscan_semaphor| %8.3f | %8.3f\n", stopc(), stop());
    // dbscan_simple(n,d,epsilon,11,array);
    // printf("dbscan_simple| %8.3f | %8.3f\n", stopc(), stop());


    //    stop(); stopc();
//    dbscan_simple(n,d,epsilon,10,array);
//    printf("dbscan simple| %8.3f | %8.3f\n", stopc(), stop());


//    for(int i=0 ; i<200 ; i++)printf("%d",processed[i]);
//    int ct=0 ;
//    for(int i=0; i<n ; i++)ct += processed[i] ;
//    printf("\ncounter=%d cores=%d\n",ct,crs);
    //joinCanonicVerbose(20, d, epsilon, array);
    // TEST TEST TEST
//    for(int i=0 ; i<n ; i++)
//        for (int j=8; j<64 ; j+=1)
//            memcpy(array+d*i+j, array+d*i, 64);
//    for(int i=0 ; i<256 ; i+=8)printf("%lf %lf %lf %lf %lf %lf %lf %lf\n", array[i],array[i+1],array[i+2],array[i+3],array[i+4],array[i+5],array[i+6],array[i+7]);
//    for(int i=n*d-256 ; i<n*d ; i+=8)printf("%lf %lf %lf %lf %lf %lf %lf %lf\n", array[i],array[i+1],array[i+2],array[i+3],array[i+4],array[i+5],array[i+6],array[i+7]);

    // dbscan(n,d,epsilon,10,array);
    // printf("dbscan| %8.3f | %8.3f\n", stopc(), stop());
    // exit(0);
    // long long * counter = corePoints(n,d,epsilon,array);
    // for (int i=0 ; i<100 ; i++) printf("%ld ", counter[i]);
    // int all=0 ; for (int i=0 ; i<n ; i++)all+=counter[i];
    // printf("\nall: %d\n",all);
    // printf("Core Points| %8.3f | %8.3f\n", stopc(), stop());
    printf("test_ego_loop");
    // test_ego_loop(n,d,epsilon,array);
    // printf("test loop 1| %8.3f | %8.3f\n", stopc(), stop());
    printf("test_ego_loop3");
    test_ego_loop3(n,d,epsilon,array);
    printf("test loop 3| %8.3f | %8.3f\n", stopc(), stop());
//    test_ego_loop2(n,d,epsilon,array);
//    printf("test loop 2| %8.3f | %8.3f\n", stopc(), stop());
//    long long * counter = corePoints(n,d,epsilon,array);
//    for (int i=0 ; i<100 ; i++) printf("%ld ", counter[i]);
//    int all=0 ; for (int i=0 ; i<n ; i++)all+=counter[i];
//    printf("\nall: %d\n",all);
//    printf("Core Points| %8.3f | %8.3f\n", stopc(), stop());

//     transClosureDens(n,d,epsilon,3,array,counter);
//     printf("TransClosur| %8.3f | %8.3f\n", stopc(), stop());
//     for(int i=0 ; i<200 ; i++)
//         printf("%5d %5d %5d\n", i, listArray[i]->minid, listArray[i]->n) ;
//     epsilonGridJoinCanonicAVX512(n, d, epsilon, array);
//     printf("EGOJoin 512| %8.3f | %8.3f\n", stopc(), stop());
//     double *self = mallocA64(sizeof (double) * n);
//     int numStripes = 20;
//     int **lower = (int **) malloc (numStripes*sizeof(int*));
//     int **upper = (int **) malloc (numStripes*sizeof(int*));
//     prepareTwoStripes(n, d, epsilon, array, lower, upper, self);
// //    for (int i=0 ; i<nn/1024 ; i++)
// //        printf("%8d %8d %8d %8d %8d\n", i, lower[0][i+nn/1024], upper[0][i+nn/1024], lower[1][i+nn/1024], upper[1][i+nn/1024]);
//     printf("Prepare    | %8.3f | %8.3f\n", stopc(), stop());
//     prepareStripes(n, d, numStripes, epsilon, array, lower, upper, self);
// //    for (int i=0 ; i<200 ; i++)
// //        printf("%8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d\n", i, lower[0][i+nn/1024], upper[0][i+nn/1024], lower[1][i+nn/1024],
// //                upper[1][i+nn/1024], lower[2][i+nn/1024], upper[2][i+nn/1024], lower[3][i+nn/1024], upper[3][i+nn/1024], lower[4][i+nn/1024], upper[4][i+nn/1024]);
// //    for (int i=0 ; i<200 ; i++)
// //        printf("%8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d\n", i, lower[0][i+nn/8], upper[0][i+nn/8], lower[1][i+nn/8],
// //                upper[1][i+nn/8], lower[2][i+nn/8], upper[2][i+nn/8], lower[3][i+nn/8], upper[3][i+nn/8], lower[4][i+nn/8], upper[4][i+nn/8]);
//     printf("PreparePar | %8.3f | %8.3f\n", stopc(), stop());
//     epsilonGridJoinSimplified(n, d, epsilon, array, numStripes, lower, upper, self);
// //    for(int i=0 ; i<20 ; i++)
// //        printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
// //                testarray[i][0], testarray[i][1], testarray[i][2], testarray[i][3], testarray[i][4], testarray[i][5], testarray[i][6], testarray[i][7], testarray[i][8], testarray[i][9],
// //                testarray[i][10], testarray[i][11], testarray[i][12], testarray[i][13], testarray[i][14], testarray[i][15], testarray[i][16], testarray[i][17], testarray[i][18], testarray[i][19]);
// //    for (int i=0 ; i<25000 ; i++){
// //        for (int s=0 ; s<numStripes ; s++)
// //            for (int j=lower[s][i+nn/8]; j<upper[s][i+nn/8] ; j++)
// //                if(testarray[i][j] != 1)
// //                    printf("ERROR ijs %d %d %d %d lower %d upper %d \n", i, j, s, testarray[i][j], lower[s][i+nn/8], upper[s][i+nn/8]) ;
// //        for (int s=1 ; s<numStripes ; s++)
// //            for (int j=upper[s-1][i+nn/8]; j<lower[s][i+nn/8] ; j++)
// //                if(testarray[i][j] != 0)
// //                    printf("ERROR ijs %d %d %d %d\n", i, j, s, testarray[i][j]) ;
// //    }
//     printf("New Simple | %8.3f | %8.3f\n", stopc(), stop());
//     JoinCanonicAVX512(n, d, epsilon, array);
//     printf("Join AVX512| %8.3f | %8.3f\n", stopc(), stop());
//     joinCanonicAVX1(n, d, epsilon, array);
//     printf("Join AVX1  | %8.3f | %8.3f\n", stopc(), stop());
//     joinCanonic(n, d, epsilon, array);
    return 0;
}
