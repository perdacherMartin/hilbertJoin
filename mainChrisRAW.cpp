/*
 * Version: Email 25. Juli
 */
 /*
  * File:   main.cpp
  * Author: boehm
  *
  * Created on 17. März 2017, 14:13
  */

 using namespace std;
 void *__gxx_personality_v0;

 #ifndef NUM_THREADS
 #define NUM_THREADS 1
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

 double timestamp(void){
     struct timeval c;
     gettimeofday(&c, NULL);
     int r = (c.tv_sec * 1000 + c.tv_usec / 1000) % 1000000 ;
     return (double)r / 1000.;
 }

 # define mytestdefine\
         _Pragma("omp parallel for")\
         for(int par=0; par<NUM_THREADS ; par++)\
             for(int i=par*n/NUM_THREADS ; i<(par+1)*n/NUM_THREADS ; i++)

 void testprag(){
     int n=1000000;
     mytestdefine{
         if(i%1000==0)printf("%d %d\n", par, i);
     }
     printf("Ready.\n");
     exit(0);
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
 typedef struct UfEl UfEl;
 struct UfEl {
     int n;
     int id;
     omp_lock_t lock;
     UfEl* parent;
 } ;
 TreeEl ** tentArray;
 ListEl ** listArray;
 UfEl * ufArray;

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

 UfEl * findUfEl(UfEl *i){
     if(i != (UfEl *)0) while (i->parent != (UfEl *) 0)
         i = i->parent;
     return i;
 }
 UfEl * findLockUfEl(UfEl *i){
     i = findUfEl(i) ;
     omp_set_lock(& i->lock) ;
     while(i->parent != (UfEl *)0){
         omp_unset_lock(& i->lock) ;
         i = findUfEl (i) ;
         omp_set_lock(& i->lock) ;
     }
     return i;
 }
 int findTestUfEl(UfEl **i, int wait){
     *i = findUfEl(*i) ;
     if(wait)
         omp_set_lock(&(*i)->lock);
     else
         if(! omp_test_lock(&(*i)->lock))
             return (0) ;
     while((*i)->parent != (UfEl *) 0){
         omp_unset_lock(&(*i)->lock) ;
         *i = findUfEl (*i) ;
         if(wait)
             omp_set_lock(&(*i)->lock);
         else
             if(! omp_test_lock(&(*i)->lock))
                 return (0) ;
     }
     return 1;
 }
 int findLockBothUfEl(UfEl **i, UfEl **j){
     for(;;) {
         findTestUfEl (i, 1) ;
         if (findTestUfEl (j, 0))
             return 0 ;
         if (*i == *j)
             return 1 ;
         omp_unset_lock (&(*i)->lock) ;
         findTestUfEl (j, 1) ;
         if (findTestUfEl (i, 0))
             return 0 ;
         if (*i == *j)
             return 1 ;
         omp_unset_lock (&(*j)->lock) ;
     }
 }
 void print_path(int i){
     printf("path %d: ", i);
     for(UfEl * h=ufArray+i ; h!= (UfEl*) 0 ; h=h->parent)
         printf("%d(%d) ", h->id, h->n);
     printf("\n");
 }
 void pathShorten(UfEl * myEh, UfEl * myEi){
  //   return;
         while (myEh != (UfEl *) 0 && myEh->parent != (UfEl *)0){
             UfEl * h = myEh ;
             myEh = myEh->parent ;
             h->parent = myEi ;
         }
 }
 void pathShorten1(UfEl * myEh, UfEl * myEi){
         while (myEh != (UfEl *) 0 && myEh->parent != (UfEl *)0){
             UfEl * h = myEh ;
             myEh = myEh->parent ;
             h->parent = myEi ;
         }
 }
 void unifyUfElDbscan(int i, int j, int iCore, int jCore) {
     // merge if both are core objects or one is core object and the other is solitaire (n==1)
     UfEl * myEi = ufArray + i;
     UfEl * myEj = ufArray + j;
     myEi = findUfEl(myEi) ;
     myEj = findUfEl(myEj);
     if (myEi == myEj) {
         if (omp_test_lock(&myEi->lock)) {
             if (myEi->parent == (UfEl*) 0) {
                 pathShorten(ufArray + i, myEi);
                 pathShorten(ufArray + j, myEj);
             }
             omp_unset_lock(&myEi->lock);
         }
         return;
     }
 //    while (myEi->parent != (UfEl *)0)
 //        myEi = myEi->parent ;
 //    while (myEj->parent != (UfEl *)0)
 //        myEj = myEj->parent ;
 //    if(myEi == myEj){
 //        if (omp_test_lock(&myEj->lock)){
 //            pathShorten(ufArray + j, myEj);
 //            pathShorten(ufArray + i, myEj);
 //            omp_unset_lock(&myEj->lock);
 //        }
 //        return;
 //    }
 //    lock1: omp_set_lock(&myEi->lock);
 //    if (myEi->parent != (UfEl *) 0){
 //        omp_unset_lock(&myEi->lock);
 //        while (myEi->parent != (UfEl *)0)
 //            myEi = myEi->parent ;
 //        goto lock1;
 //    }
 //    pathShorten(ufArray + i, myEi);
 //    lock2:
 //    if(myEi == myEj){
 //        pathShorten(ufArray + j, myEj);
 //        omp_unset_lock(&myEi->lock);
 //        return;
 //    }
 //    if (!omp_test_lock(&myEj->lock)){
 //        omp_unset_lock(&myEi->lock);
 //        goto lock3 ;
 //    } else {
 //        if(myEj->parent != (UfEl *)0){
 //            omp_unset_lock(&myEj->lock);
 //            while (myEj->parent != (UfEl *)0)
 //                myEj = myEj->parent ;
 //            goto lock2 ;
 //        }
 //        goto locksOk;
 //    }
 //    lock3:
 //    omp_set_lock(&myEj->lock);
 //    if (myEj->parent != (UfEl *) 0){
 //        omp_unset_lock(&myEj->lock);
 //        while (myEj->parent != (UfEl *)0)
 //            myEj = myEj->parent ;
 //        goto lock3;
 //    }
 //    pathShorten(ufArray + j, myEj);
 //    lock4:
 //    if(myEi == myEj){
 //        pathShorten(ufArray + i, myEi);
 //        omp_unset_lock(&myEj->lock);
 //        return;
 //    }
 //    if (!omp_test_lock(&myEi->lock)){
 //        omp_unset_lock(&myEj->lock);
 //        goto lock1 ;
 //    } else {
 //        if(myEi->parent != (UfEl *)0){
 //            omp_unset_lock(&myEi->lock);
 //            while (myEi->parent != (UfEl *)0)
 //                myEi = myEi->parent ;
 //            goto lock4 ;
 //        }
 //            // goto locksOk;
 //    }
     if(findLockBothUfEl(&myEi, &myEj)){
         pathShorten(ufArray+i, myEi);
         pathShorten(ufArray+j, myEi);
         omp_unset_lock(&myEi->lock);
         return;
     }
     locksOk:
     if (iCore && (jCore || myEj->n == 1) || jCore && myEi->n == 1) {
         if (myEi->n < myEj->n) {
             myEi->parent = myEj ;
             myEj->n += myEi->n ;
             pathShorten(ufArray+i, myEj) ;
             pathShorten(ufArray+j, myEj) ;
         } else {
             myEj->parent = myEi ;
             myEi->n += myEj->n ;
             pathShorten(ufArray+i, myEi) ;
             pathShorten(ufArray+j, myEi) ;
         }
     }
     //    printf("%2d UNLOCK %2d %2d\n", count, mylock1->minid, mylock2->minid);
     omp_unset_lock(&myEi->lock);
     omp_unset_lock(&myEj->lock);
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
             loadstart[i] = loadstart[i - 1];\
     loadstart[NUM_THREADS+1] = loadstart[NUM_THREADS];

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
 #define KBLOCK 8
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

 int *reorder_dim;

 void reorder_dimensions(int n, int d, double *array) {
     if (d > 128) {
 #pragma omp parallel for
         for (int par = 0; par < NUM_THREADS; par++) {
             double help[d + 8]__attribute__((aligned(64)));
             for (int i = par * n / NUM_THREADS; i < (par + 1) * n / NUM_THREADS; i++) {
                 memcpy(help, array + i*d, d * sizeof (double));
                 for (int j = 0; j < d; j++)
                     array[reorder_dim[j]+i*d] = help[j];
             }
         }
         return;
     }
 #pragma omp parallel for
     for (int par = 0; par < NUM_THREADS; par++) {
         __m256i i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15;
         __m256i * rd = (__m256i *) reorder_dim;
         switch ((d - 1) / 8) {
             case 15: i15 = _mm256_loadu_si256(rd++); //rd += 8;
             case 14: i14 = _mm256_loadu_si256(rd++); //rd += 8;
             case 13: i13 = _mm256_loadu_si256(rd++); //rd += 8;
             case 12: i12 = _mm256_loadu_si256(rd++); //rd += 8;
             case 11: i11 = _mm256_loadu_si256(rd++); //rd += 8;
             case 10: i10 = _mm256_loadu_si256(rd++); //rd += 8;
             case 9: i9 = _mm256_loadu_si256(rd++); //rd += 8;
             case 8: i8 = _mm256_loadu_si256(rd++); //rd += 8;
             case 7: i7 = _mm256_loadu_si256(rd++); //rd += 8;
             case 6: i6 = _mm256_loadu_si256(rd++); //rd += 8;
             case 5: i5 = _mm256_loadu_si256(rd++); //rd += 8;
             case 4: i4 = _mm256_loadu_si256(rd++); //rd += 8;
             case 3: i3 = _mm256_loadu_si256(rd++); //rd += 8;
             case 2: i2 = _mm256_loadu_si256(rd++); //rd += 8;
             case 1: i1 = _mm256_loadu_si256(rd++); //rd += 8;
             case 0: i0 = _mm256_loadu_si256(rd);
         }
         __mmask8 k = 256 - (1 << ((128 - d) % 8));
         vec x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15;
         for (int i = par * n / NUM_THREADS; i < (par + 1) * n / NUM_THREADS; i++) {
             double * xd = array + i*d;
             switch ((d - 1) / 8) {
                 case 15: x15 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 14: x14 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 13: x13 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 12: x12 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 11: x11 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 10: x10 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 9: x9 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 8: x8 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 7: x7 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 6: x6 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 5: x5 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 4: x4 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 3: x3 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 2: x2 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 1: x1 = _mm512_loadu_pd(xd);
                     xd += 8;
                 case 0: x0 = _mm512_loadu_pd(xd);
             }
             xd = array + i*d;
             switch ((d - 1) / 8) {
                 case 15: _mm512_i32scatter_pd(xd, i15, x15, 8);
                 case 14: _mm512_i32scatter_pd(xd, i14, x14, 8);
                 case 13: _mm512_i32scatter_pd(xd, i13, x13, 8);
                 case 12: _mm512_i32scatter_pd(xd, i12, x12, 8);
                 case 11: _mm512_i32scatter_pd(xd, i11, x11, 8);
                 case 10: _mm512_i32scatter_pd(xd, i10, x10, 8);
                 case 9: _mm512_i32scatter_pd(xd, i9, x9, 8);
                 case 8: _mm512_i32scatter_pd(xd, i8, x8, 8);
                 case 7: _mm512_i32scatter_pd(xd, i7, x7, 8);
                 case 6: _mm512_i32scatter_pd(xd, i6, x6, 8);
                 case 5: _mm512_i32scatter_pd(xd, i5, x5, 8);
                 case 4: _mm512_i32scatter_pd(xd, i4, x4, 8);
                 case 3: _mm512_i32scatter_pd(xd, i3, x3, 8);
                 case 2: _mm512_i32scatter_pd(xd, i2, x2, 8);
                 case 1: _mm512_i32scatter_pd(xd, i1, x1, 8);
                 case 0: _mm512_mask_i32scatter_pd(xd, k, i0, x0, 8);
             }
         }
     }
 }

 static inline void transpose_8xd_reorder(int n, int d, double *EGO_array) {
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
                 _mm512_store_pd(help + reorder_dim[k] * 8, v1);
                 _mm512_store_pd(help + reorder_dim[k+1] * 8, v2);
                 _mm512_store_pd(help + reorder_dim[k+2] * 8, v3);
                 _mm512_store_pd(help + reorder_dim[k+3] * 8, v4);
                 _mm512_store_pd(help + reorder_dim[k+4] * 8, v5);
                 _mm512_store_pd(help + reorder_dim[k+5] * 8, v6);
                 _mm512_store_pd(help + reorder_dim[k+6] * 8, v7);
                 _mm512_store_pd(help + reorder_dim[k+7] * 8, v8);
             }
             memcpy(EGO_array + i * 8 * d, help, d * 8 * sizeof (double));
         }
     }
 }

 static inline void transpose_dx8_reorder(int n, int d, double *EGO_array) {
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
                 vec v1 = _mm512_load_pd(EGO_array + 8 * (i*d + reorder_dim[k]));
                 vec v2 = _mm512_load_pd(EGO_array + 8 * (i*d + reorder_dim[k+1]));
                 vec v3 = _mm512_load_pd(EGO_array + 8 * (i*d + reorder_dim[k+2]));
                 vec v4 = _mm512_load_pd(EGO_array + 8 * (i*d + reorder_dim[k+3]));
                 vec v5 = _mm512_load_pd(EGO_array + 8 * (i*d + reorder_dim[k+4]));
                 vec v6 = _mm512_load_pd(EGO_array + 8 * (i*d + reorder_dim[k+5]));
                 vec v7 = _mm512_load_pd(EGO_array + 8 * (i*d + reorder_dim[k+6]));
                 vec v8 = _mm512_load_pd(EGO_array + 8 * (i*d + reorder_dim[k+7]));
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

 //    omp_lock_t loadlocks[EGO_stripes]; for(int s=0; s<EGO_stripes; s++) omp_init_lock(loadlocks+s);\

 # define EGO_PARALLEL_TRAN(n,d,epsilon,stripes,array) {\
     int EGO_n = (n);\
     int EGO_d = (d);\
     double EGO_epsilon = (epsilon);\
     double * EGO_array = (array);\
     int EGO_blocks = (EGO_d + KBLOCK - 1) / KBLOCK;\
     int EGO_stripes = (stripes);\
     unsigned long long usedload[NUM_THREADS]; for(int i=0; i<NUM_THREADS; i++) usedload[i]=0ull;\
     omp_lock_t criticallock; omp_init_lock(&criticallock); int scritical = -1;\
     epsilonGridOrdering(EGO_n, EGO_d, EGO_epsilon, EGO_array);\
     int nn = ceilpowtwo(EGO_n);\
     int **lower = (int **) malloc (EGO_stripes*sizeof(int*));\
     int **upper = (int **) malloc (EGO_stripes*sizeof(int*));\
     double *self = callocA64(sizeof (double) * EGO_n * EGO_blocks);\
     prepareStripes(EGO_n, EGO_d, EGO_stripes, EGO_epsilon, EGO_array, lower, upper, (double *)0);\
     EGO_epsilon = EGO_epsilon * EGO_epsilon / 2;\
     for(int i=EGO_n ; i<(EGO_n+7)/8*8 ; i++)\
         for(int j=0 ; j<EGO_d ; j++)\
             EGO_array[i*EGO_d+j] = 1e150 - 1e140*(double)(i%8);\
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
     unsigned long long *sloadcumcum = (unsigned long long *) malloc(EGO_stripes * (n+7)/8 * sizeof(unsigned long long));\
     sloadcumcum[EGO_stripes-1]=upper[EGO_stripes-1][nn/8]-lower[EGO_stripes-1][nn/8];\
     for(int s=EGO_stripes-1; s>=0 ; s--)\
         sloadcumcum[s] = sloadcumcum[s+1] + upper[s][nn/8] - lower[s][nn/8];\
     for(int i=1; i<(n+7)/8 ; i++)\
         for(int s=EGO_stripes-1; s>=0 ; s--)\
             sloadcumcum[i*EGO_stripes+s] = sloadcumcum[(i-1)*EGO_stripes+s] + upper[s][nn/8+i] - lower[s][nn/8+i];\
     for(int i=0; i<(n+7)/8 ; i++)\
         for(int s=EGO_stripes-2; s>=0 ; s--)\
             sloadcumcum[i*EGO_stripes+s] += sloadcumcum[i*EGO_stripes+s+1];\
     unsigned long long * assignedload = (unsigned long long *) calloc(NUM_THREADS, sizeof (unsigned long long));\
     long long overall_load = 0;\
     for (int i = 0; i < EGO_n / 8; i++)\
         for (int j = 0; j < EGO_stripes; j++)\
             overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
     int *loadstart = (int *) calloc((NUM_THREADS + 1) * EGO_stripes, sizeof(int));\
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

 //            if(omp_test_lock(loadlocks+s)) {\
 //                omp_set_lock(&criticallock);\
 //                printf("Thread %d (par=%d) saving the load %ld\n", omp_get_thread_num(), par, usedload[s]);\
 //                memcpy(savedload+NUM_THREADS*s, usedload, NUM_THREADS*sizeof(long long)) ;\
 //                scritical = s;\
 //                omp_unset_lock(&criticallock);\
 //            } else while (s>scritical){omp_set_lock(&criticallock); omp_unset_lock(&criticallock);}\

 int cmp_ull(const void *a, const void *b){
     if(*(unsigned long long *)a  >  *(unsigned long long *)b) return 1;
     if(*(unsigned long long *)a  <  *(unsigned long long *)b) return -1;
     return 0;
 }
 inline double limitLoad(int m, unsigned long long *row, unsigned long long remainder){
     unsigned long long h=(remainder+m-1)/m;
     for(int i=0; i<m ; i++)
         if(row[i]>h) goto strongtest;
     h=remainder;
     for(int i=0; i<m ; i++)
         h+=row[i];
     return(double)h/(double)m ;
 strongtest:
     unsigned long long rowcp[m];
     memcpy(rowcp, row, m*sizeof(unsigned long long));
     qsort(rowcp, m, sizeof(unsigned long long), cmp_ull);
 //    printf("%ld %ld %ld %ld\n", rowcp[0], rowcp[1], rowcp[2], rowcp[3]);
     for(int i=0; i<m-1 ; i++){
         remainder += rowcp[i] ;
         if (rowcp[i+1] >= (remainder+i)/(i+1))
             return (double)remainder/(i+1);
     }
     return (double)(remainder+rowcp[m-1])/m;
 }

 #define EGO_LOOP_TRAN_LOADBAL\
         for (int s = 0; s < EGO_stripes; s++) {\
             if(s > scritical - 1){\
                 if(!omp_test_lock(&criticallock)){printf ("WAIT %d\n", par);omp_set_lock(&criticallock);}\
                 if(s > scritical){\
                     printf("Thread %d (par=%d) saving the load %ld\n", omp_get_thread_num(), par, usedload[par]);\
                     unsigned long long hangover[NUM_THREADS];\
                     memcpy(hangover, usedload, NUM_THREADS*sizeof(long long)) ;\
                     for(int i=0; i<NUM_THREADS ; i++)\
                         hangover[i] = assignedload[i]-hangover[i];\
                     double stdLoad = limitLoad(NUM_THREADS, hangover, sloadcumcum[(n-1)/8*EGO_stripes+s]);\
                     printf("stdLoad = %f\n", stdLoad);\
                     unsigned long long cumHangover = 0ull;\
                     int curN = 0;\
                     for(int i=0; i<NUM_THREADS ; i++){\
                         for(int stepsize=n/2/NUM_THREADS ; stepsize >= 4 ; stepsize /= 4)\
                             while(curN+stepsize<(n+7)/8 && (double) (sloadcumcum[(curN+stepsize)*EGO_stripes+s]) < (double)i*stdLoad-(double)cumHangover) curN+=stepsize;\
                         while(curN<(n+7)/8 && (double) (sloadcumcum[curN*EGO_stripes+s]) < (double)i*stdLoad-(double)cumHangover) curN++;\
                         cumHangover += min(hangover[i],stdLoad) ;\
                         /*printf("%d (%d %ld)  ", curN, loadstart[s*(NUM_THREADS+1)+i], cumHangover);*/\
                         loadstart[s*(NUM_THREADS+1)+i]=curN ;\
                         if(i>0 && curN>0) {unsigned long long h=sloadcumcum[s+(curN-1)*EGO_stripes]-(s+1<EGO_stripes?sloadcumcum[s+1+(curN-1)*EGO_stripes]:0); assignedload[i-1]+=h; assignedload[i]-=h;}\
                     }\
                     loadstart[s*(NUM_THREADS+1)+NUM_THREADS]=(n+7)/8;\
                     assignedload[NUM_THREADS-1]+=sloadcumcum[s+(n-1)/8*EGO_stripes]-(s+1<EGO_stripes?sloadcumcum[s+1+(n-1)/8*EGO_stripes]:0);\
                     /*printf("%d\n", (n+7)/8);*/\
                     /*printf("ass %ld %ld %ld %ld %ld %ld %ld %ld\n",assignedload[0],assignedload[1],assignedload[2],assignedload[3],assignedload[4],assignedload[5],assignedload[6],assignedload[7]);*/\
                     scritical = s;\
                 }\
                 omp_unset_lock(&criticallock);\
             }\
             int imin = loadstart[s*(NUM_THREADS+1)+par];\
             int imax = loadstart[s*(NUM_THREADS+1)+par+1];\
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
                 usedload[par]++;\
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
     long long result = 0;
     long long refinements = 0;
     unsigned long long savedload[5*NUM_THREADS];
     double starttimestamp = timestamp() ;
     EGO_PARALLEL_TRAN(n, d, epsilon, 41, array)
         printf("timestamp index ready %6.2f\n",timestamp()-starttimestamp);
 //        for(int i=0 ; i<NUM_THREADS + 4 ; i+=4)
 //            printf("%2d %9d %9d %9d %9d\n", i, loadstart[i], loadstart[i+1], loadstart[i+2], loadstart[i+3]);
 //        for (int x=0 ; x<NUM_THREADS ; x++){
 //            long long cum_load = 0;
 //            for(int i = loadstart[x] ; i<loadstart[x+1] ; i++)
 //                for (int j = 0; j < EGO_stripes; j++)
 //                    cum_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];
 //            if (x%4 == 0) printf("\n %2d", x);
 //            printf("%9ld (%5.2f%%)  ", cum_load, 100.*(double)cum_load/(double)overall_load);
 //        }
         printf("overall_load: %ld / %ld (=n*(n-1)/2 / 64) ==> %f\n", overall_load, (long long)n/128*(n-1), (double)overall_load/n/(n-1)*128);

        #pragma omp parallel for proc_bind(close) reduction(+:result) reduction(+:refinements)
     EGO_PREPARE
         veci64 resultvec = _mm512_setzero_si512();
         veci64 eights = _mm512_set1_epi64(8ll) ;
         long long refineload = 0;
     EGO_LOOP_TRAN {
         resultvec += eights - _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63) ;
         refineload ++;
     }
     EGO_CONSOLIDATE{
         result += _mm512_reduce_add_epi64(resultvec);
         refinements += refineload ;
         int curload=0;
         for(int i=loadstart[par] ; i<loadstart[par+1] ; i++)
             for(int s=0 ; s<EGO_stripes ; s++)
                 curload += upper[s][i+nn/8] - lower[s][i+nn/8];
         printf("Consolidate %6.2f %d %d %d %d %d %ld %ld\n",timestamp()-starttimestamp, par, omp_get_thread_num(), loadstart[par], loadstart[par+1]-loadstart[par], curload, refineload, result);

 //        double testres[8] __attribute__((aligned(64)));
 //        _mm512_store_epi64(testres, resultvec);
 //        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
     }
     EGO_END_TRAN
     printf("result %ld\n", result);
     printf("refinements %ld Mio (%ld 8x8)\n", refinements*64/1000000, refinements);
 //    for(int par=0 ; par<NUM_THREADS ; par++, printf("\n"))
 //        for(int s=0 ; s<5 ; s++)
 //            printf("%ld ",savedload[NUM_THREADS*s+par]);

 }

 int test_ego_loop_nested(int n, int d, double epsilon, double *array){
     long long result = 0;
     unsigned long long savedload[5*NUM_THREADS];
     double starttimestamp = timestamp() ;
     omp_set_nested(1);
     EGO_PARALLEL_TRAN(n, d, epsilon, 5, array)
         #pragma omp parallel for num_threads(4) proc_bind(spread)
 //    EGO_PREPARE
 for (int p1=0 ; p1<4 ; p1++)
 #pragma omp parallel for num_threads(16) proc_bind(close)
 for (int par=p1*16; par<(p1+1)*16 ; par++){
         veci64 resultvec = _mm512_setzero_si512();
         veci64 eights = _mm512_set1_epi64(8ll) ;
         long long refineload = 0;
     EGO_LOOP_TRAN {
         resultvec += eights - _mm512_srli_epi64(_mm512_castpd_si512(sum1), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum2), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum3), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum4), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum5), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum6), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum7), 63)
                             - _mm512_srli_epi64(_mm512_castpd_si512(sum8), 63) ;
         refineload ++;
     }
     EGO_CONSOLIDATE{
 #pragma omp critical
         {
         result += _mm512_reduce_add_epi64(resultvec);
         }
         int curload=0;
         for(int i=loadstart[par] ; i<loadstart[par+1] ; i++)
             for(int s=0 ; s<EGO_stripes ; s++)
                 curload += upper[s][i+nn/8] - lower[s][i+nn/8];
         printf("Consolidate %6.2f %d %d %d %d %d %d %ld %ld\n",timestamp()-starttimestamp, p1, omp_get_thread_num(), par, loadstart[par], loadstart[par+1]-loadstart[par], curload, refineload, result);

 //        double testres[8] __attribute__((aligned(64)));
 //        _mm512_store_epi64(testres, resultvec);
 //        printf("par = %d: %d %d\n", par, result, testres[0]+testres[1]+testres[2]+testres[3]+testres[4]+testres[5]+testres[6]+testres[7]);
     }
     EGO_END_TRAN
     printf("result %ld\n", result);
     for(int par=0 ; par<NUM_THREADS ; par++, printf("\n"))
         for(int s=0 ; s<5 ; s++)
             printf("%ld ",savedload[NUM_THREADS*s+par]);

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
     EGO_PARALLEL_TRAN(n, d, epsilon, 5, array)
         #pragma omp parallel for
     EGO_PREPARE
         long long * locct = (long long*)callocA64(n*sizeof(long long));
 //        long long * temp = (long long*)mallocA64(64*sizeof(long long));
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
         printf("Consolidate thread %d\n",omp_get_thread_num());
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
 long long * corePointsMutex(int n, int d, double epsilon, double *array){
     long long * counter = (long long*)callocA64(n*sizeof(long long));
     omp_lock_t locks[(n+7)/8] ;
     for (int i=0 ; i<(n+7)/8 ; i++)
         omp_init_lock(locks+i);
     // Aufruf EGO-Join EGO_PARALLEL_TRAN ... EGO_END_TRAN
     EGO_PARALLEL_TRAN(n, d, epsilon, 2, array)
         #pragma omp parallel for
     EGO_PREPARE
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
         omp_set_lock(locks+i);
         _mm512_store_si512(counter + i*8, _mm512_load_si512(counter + i*8) + eights - (vi1+vi2+vi3+vi4+vi5+vi6+vi7+vi8)) ;
         omp_unset_lock(locks+i);
         transposeAVX512i(vi1,vi2,vi3,vi4,vi5,vi6,vi7,vi8);
         omp_set_lock(locks+j);
         _mm512_store_si512(counter + j*8, _mm512_load_si512(counter + j*8) + eights - (vi1+vi2+vi3+vi4+vi5+vi6+vi7+vi8)) ;
         omp_unset_lock(locks+j);
     }
     EGO_CONSOLIDATE
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
     double starttimestamp=timestamp();
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
         printf("Consolidate point count %f %d (thread=%d)\n",timestamp()-starttimestamp, par, omp_get_thread_num());
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
         printf("Consolidate thread %f %d\n",timestamp()-starttimestamp,par);
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
 #define BLA(x,y)  \
         _mm512_store_pd(temp, x);\
         for (int k = 0; k < 8; k++)\
             if (temp[k] >= 0) {\
                 if (locct[par][j * 8 + y] < minPts) {\
                     memory[par][(j * 8 + y) * minPts + locct[par][j * 8 + y]++] = i * 8 + k;\
                     if (locct[par][i * 8 + k] < minPts)\
                         memory[par][(i * 8 + k) * minPts + locct[par][i * 8 + k]] = j * 8 + y;\
                     locct[par][i * 8 + k]++;\
                 } else {\
                     if (locct[par][i * 8 + k] < minPts)\
                         memory[par][(i * 8 + k) * minPts + locct[par][i * 8 + k]++] = j * 8 + y;\
                     else{\
                         unifyList(i * 8 + k, j * 8 + y, 0);num_unify++;\
                         locct[par][i * 8 + k]++;locct[par][j * 8 + y]++;\
             }   }   }

 //load_unify += min(listArray[i*8+k]->n, listArray[j*8+y]->n);
 //num_unify++;
 void dbscanOnepass(int n, int d, double epsilon, int minPts, double *array){
     int stripes = 5;
     double starttimestamp=timestamp();
     unsigned long long *savedload = (unsigned long long *) malloc(NUM_THREADS*stripes*sizeof(unsigned long long));
     double eps_sq = epsilon * epsilon ;
     long long * counter = (long long*)callocA64(n*sizeof(long long));
     long long * locct[NUM_THREADS] ;
     int * memory[NUM_THREADS] ;
     omp_lock_t locks[2*NUM_THREADS] ;
     for (int i=0 ; i<2*NUM_THREADS ; i++)
         omp_init_lock(locks+i);
     listArray = (ListEl **) malloc(n * sizeof (ListEl*));
     listArray[0] = (ListEl *) malloc(n * sizeof (ListEl));
     for (int i = 0; i < n; i++) {
         listArray[i] = listArray[0] + i;
         listArray[i]->n = 1;
         listArray[i]->minid = i;
         listArray[i]->left = (ListEl*) 0;
         omp_init_lock(&listArray[i]->lock);
     }
     EGO_PARALLEL_TRAN(n, d, epsilon, stripes, array)
         #pragma omp parallel for
     EGO_PREPARE
         long long num_unify = 0;
         long long load_unify = 0;
         locct[par] = (long long*)callocA64(n*sizeof(long long));
         memory[par] = (int *)mallocA64(n*minPts*sizeof(int));
         double temp[8]__attribute__((aligned(64)));
     EGO_LOOP_TRAN_LOADBAL{
         load_unify++;
         BLA(sum1, 0);
         BLA(sum2, 1);
         BLA(sum3, 2);
         BLA(sum4, 3);
         BLA(sum5, 4);
         BLA(sum6, 5);
         BLA(sum7, 6);
         BLA(sum8, 7);
     }
     EGO_CONSOLIDATE{
         int curload=0;
         for(int i=loadstart[par] ; i<loadstart[par+1] ; i++)
             for(int s=0 ; s<EGO_stripes ; s++)
                 curload += upper[s][i+nn/8] - lower[s][i+nn/8];
         printf("Consolidate %6.2f %d %d %d %d %ld %ld\n",timestamp()-starttimestamp, par, loadstart[par], loadstart[par+1]-loadstart[par], curload, num_unify, load_unify);
         int segments_added=0;
         int curseg = 2*(omp_get_thread_num()%NUM_THREADS) ;
         char * segment=(char*) calloc(2*NUM_THREADS,1);
         while(segments_added < 2 * NUM_THREADS){
             if(!segment[curseg] && omp_test_lock(locks+curseg)){
                 for(int i=curseg*n/NUM_THREADS/16*8 ; i<(curseg+1)*n/NUM_THREADS/16*8 ; i+=8)
                     _mm512_store_si512(counter+i, _mm512_add_epi64(_mm512_load_si512(counter+i),_mm512_load_si512(locct[par]+i)));
                 omp_unset_lock(locks+curseg);
                 segments_added++;
                 segment[curseg] = 1;
             }
             if(++curseg >= 2*NUM_THREADS)
                 curseg = 0;
         }
     }
     EGO_END_TRAN
     printf("Join ready. Cleanup.\n");
 #pragma omp parallel for
     for(int par=0 ; par<NUM_THREADS ; par++){
         for(int i=0 ; i<n ; i++)
             for(int j=0 ; j<minPts && j<locct[par][i] ; j++)
                 unifyListDbscan(i, memory[par][i*minPts+j], counter[i]>=minPts, counter[memory[par][i*minPts+j]]>=minPts);
     }
     for (int i = 0; i < 20; i++)printf("%ld ", counter[i]);
     printf("\n");
     for (int i = n - 20; i < n; i++)printf("%ld ", counter[i]);
     printf("\n");

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
                 if(listArray[j]->minid == i && counter[j] >= minPts) ct++ ;
             printf("%d %d(core %d) %d ", listArray[i]->minid, listArray[i]->n, ct, counter[i]) ;
             printf("(%d %d %d ...)\n", listArray[i]->left->minid, listArray[i]->left->left->minid, listArray[i]->left->left->left->minid);
         }
     for(int par=0 ; par<NUM_THREADS ; printf("%d\n", par++))
         for(int s=0 ; s<5 ; s++)
             printf("%ld ",savedload[NUM_THREADS*s+par]);
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

 #define BLAUf(x,y)  \
         _mm512_store_pd(temp, x);\
         for (int k = 0; k < 8; k++)\
             if (temp[k] >= 0) {\
                 if (locct[par][j * 8 + y] < minPts) {\
                     memory[par][(j * 8 + y) * minPts + locct[par][j * 8 + y]++] = i * 8 + k;\
                     if (locct[par][i * 8 + k] < minPts)\
                         memory[par][(i * 8 + k) * minPts + locct[par][i * 8 + k]] = j * 8 + y;\
                     locct[par][i * 8 + k]++;\
                 } else {\
                     if (locct[par][i * 8 + k] < minPts)\
                         memory[par][(i * 8 + k) * minPts + locct[par][i * 8 + k]++] = j * 8 + y;\
                     else{\
                         unifyUfElDbscan(i * 8 + k, j * 8 + y, 1,1);num_unify++;\
                         locct[par][i * 8 + k]++;locct[par][j * 8 + y]++;\
             }   }   }

 //load_unify += min(listArray[i*8+k]->n, listArray[j*8+y]->n);
 //num_unify++;
 void dbscanOnepassUf(int n, int d, double epsilon, int minPts, double *array){
     int stripes = 5;
     double starttimestamp=timestamp();
     unsigned long long *savedload = (unsigned long long *) malloc(NUM_THREADS*stripes*sizeof(unsigned long long));
     double eps_sq = epsilon * epsilon ;
     long long * counter = (long long*)callocA64(n*sizeof(long long));
     long long * locct[NUM_THREADS] ;
     int * memory[NUM_THREADS] ;
     omp_lock_t locks[2*NUM_THREADS] ;
     for (int i=0 ; i<2*NUM_THREADS ; i++)
         omp_init_lock(locks+i);
     ufArray = (UfEl *) malloc(n * sizeof (UfEl));
     for (int i = 0; i < n; i++) {
         ufArray[i].n = 1;
         ufArray[i].id = i;
         ufArray[i].parent = (UfEl*) 0;
         omp_init_lock(&ufArray[i].lock);
     }
     EGO_PARALLEL_TRAN(n, d, epsilon, stripes, array)
         #pragma omp parallel for
     EGO_PREPARE
         long long num_unify = 0;
         long long load_unify = 0;
         locct[par] = (long long*)callocA64(n*sizeof(long long));
         memory[par] = (int *)mallocA64(n*minPts*sizeof(int));
         double temp[8]__attribute__((aligned(64)));
     EGO_LOOP_TRAN_LOADBAL{
         load_unify++;
         BLAUf(sum1, 0);
         BLAUf(sum2, 1);
         BLAUf(sum3, 2);
         BLAUf(sum4, 3);
         BLAUf(sum5, 4);
         BLAUf(sum6, 5);
         BLAUf(sum7, 6);
         BLAUf(sum8, 7);
     }
     EGO_CONSOLIDATE{
         int curload=0;
         for(int i=loadstart[par] ; i<loadstart[par+1] ; i++)
             for(int s=0 ; s<EGO_stripes ; s++)
                 curload += upper[s][i+nn/8] - lower[s][i+nn/8];
         printf("Consolidate %6.2f %d %d %d %d %ld %ld\n",timestamp()-starttimestamp, par, loadstart[par], loadstart[par+1]-loadstart[par], curload, num_unify, load_unify);
         int segments_added=0;
         int curseg = 2*(omp_get_thread_num()%NUM_THREADS) ;
         char * segment=(char*) calloc(2*NUM_THREADS,1);
         while(segments_added < 2 * NUM_THREADS){
             if(!segment[curseg] && omp_test_lock(locks+curseg)){
                 for(int i=curseg*n/NUM_THREADS/16*8 ; i<(curseg+1)*n/NUM_THREADS/16*8 ; i+=8)
                     _mm512_store_si512(counter+i, _mm512_add_epi64(_mm512_load_si512(counter+i),_mm512_load_si512(locct[par]+i)));
                 omp_unset_lock(locks+curseg);
                 segments_added++;
                 segment[curseg] = 1;
             }
             if(++curseg >= 2*NUM_THREADS)
                 curseg = 0;
         }
     }
     EGO_END_TRAN
     printf("Join ready. Cleanup.\n");
     for(int iI=0 ; iI<100 ; iI++)
         print_path(iI*n/100);
     for (int iI=0 ; iI<10 ; iI++)
         print_path(iI);
     pathShorten1(ufArray, findUfEl(ufArray));
     for (int iI=0 ; iI<10 ; iI++)
         print_path(iI);
 #pragma omp parallel for
     for(int par=0 ; par<NUM_THREADS ; par++){
         for(int i=0 ; i<n ; i++)
             for(int j=0 ; j<minPts && j<locct[par][i] ; j++)
                 unifyUfElDbscan(i, memory[par][i*minPts+j], counter[i]>=minPts, counter[memory[par][i*minPts+j]]>=minPts);
     }
     for (int i = 0; i < 20; i++)printf("%ld ", counter[i]);
     printf("\n");
     for (int i = n - 20; i < n; i++)printf("%ld ", counter[i]);
     printf("\n");

     //    printf("UnifyCounter = %ld\n", unifyCounter) ;
 //    for(int i=0 ; i<HISTO_SIZE ; i++)
 //        printf("%d ", histo[i]);
 //    printf("\n");
 //    for(int i=0 ; i<20 ; i++)
 //        printf("%d %d %d %d\n", i, listArray[i]->minid, listArray[i]->n, counter[i]) ;
 //    printf("\n");
 //    for(int i=0 ; i<n ; i++)
 //        if(ufArray[i].id == i && ufArray[i].n > 1000){
 //            int ct = 0;
 //            for (int j=i ; j<n ; j++)
 //                if(ufArray[j].id == i && counter[j] >= minPts) ct++ ;
 //            printf("%d %d(core %d) %d ", ufArray[i].id, ufArray[i].n, ct, counter[i]) ;
 //            printf("(%d %d %d ...)\n", ufArray[i].parent->id, ufArray[i].parent->parent->id, ufArray[i].parent->parent->parent->id);
 //        }
     for(int par=0 ; par<NUM_THREADS ; printf("%d\n", par++))
         for(int s=0 ; s<5 ; s++)
             printf("%ld ",savedload[NUM_THREADS*s+par]);
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

 inline void submitUnifyList(int * s, int ** commField, int * commCounter, omp_lock_t testLock){
 #pragma omp critical(commCritical)
     {
         if(*commCounter == 0)
             omp_unset_lock(&testLock);
         commField[(*commCounter)++] = s;
     }
 }
 #define BLA2(x,y)  \
         _mm512_store_pd(temp, x);\
         for (int k = 0; k < 8; k++)\
             if (temp[k] >= 0) {\
                 if (locct[par][j * 8 + y] < minPts) {\
                     memory[par][(j * 8 + y) * minPts + locct[par][j * 8 + y]++] = i * 8 + k;\
                     if (locct[par][i * 8 + k] < minPts)\
                         memory[par][(i * 8 + k) * minPts + locct[par][i * 8 + k]] = j * 8 + y;\
                     locct[par][i * 8 + k]++;\
                 } else {\
                     if (locct[par][i * 8 + k] < minPts)\
                         memory[par][(i * 8 + k) * minPts + locct[par][i * 8 + k]++] = j * 8 + y;\
                     else{\
                         if(commCount>=commSize){\
                             commArea[2*commCount] = -1 ;\
                             submitUnifyList(commArea, commField, &commCounter, testlock);\
                             commSize = 1024*1024<<(commCounter/64);\
                             commArea = (int *) malloc((2*commSize+2)*sizeof(int));\
                             if(! commArea){printf("Out of memory.\n"); exit(-1);}\
                             commCount = 0;\
                         }\
                         commArea[2*commCount] = i * 8 + k;\
                         commArea[2*commCount+1] = j * 8 + y;\
                         commCount++;\
                         num_unify++;\
                         locct[par][i * 8 + k]++;locct[par][j * 8 + y]++;\
             }   }   }

 //load_unify += min(listArray[i*8+k]->n, listArray[j*8+y]->n);
 //num_unify++;
 void dbscanOnepass2(int n, int d, double epsilon, int minPts, double *array){
     omp_lock_t testlock;
     int **commField=(int **)malloc(1024*sizeof(int*));
     int commCounter = 0;
     int done = 0;
     omp_init_lock(&testlock);
 //    omp_set_nested(1);
 //    omp_set_dynamic(1);
     omp_set_lock(&testlock);
     printf("Ready.\n");
     double starttimestamp=timestamp();
     double eps_sq = epsilon * epsilon ;
     long long * counter = (long long*)callocA64(n*sizeof(long long));
     long long * locct[NUM_THREADS] ;
     int * memory[NUM_THREADS] ;
     omp_lock_t locks[2*NUM_THREADS] ;
     for (int i=0 ; i<2*NUM_THREADS ; i++)
         omp_init_lock(locks+i);
     listArray = (ListEl **) malloc(n * sizeof (ListEl*));
     listArray[0] = (ListEl *) malloc(n * sizeof (ListEl));
     for (int i = 0; i < n; i++) {
         listArray[i] = listArray[0] + i;
         listArray[i]->n = 1;
         listArray[i]->minid = i;
         listArray[i]->left = (ListEl*) 0;
         omp_init_lock(&listArray[i]->lock);
     }
     EGO_PARALLEL_TRAN(n, d, epsilon, 5, array)
     omp_set_num_threads(NUM_THREADS + 1);
 #pragma omp parallel
 #pragma omp master
     {
     for (int par = 0; par < NUM_THREADS; par++) {
 #pragma omp task
     {
         if(par==0){
 #pragma omp task
             {
                 while(done<NUM_THREADS){
                     while(commCounter==0){
                         printf("Waiting for unifications %d\n", omp_get_thread_num());
                         omp_set_lock(&testlock);
                         printf("Continuing doing unifications\n");
                     }
                     int * commArea;
 #pragma omp critical(commCritical)
                     {
                         commArea = commField[--commCounter];
                     }
                     //printf("unifying [%d %d][%d %d][%d %d] ... ", commArea[0],commArea[1],commArea[2],commArea[3],commArea[4],commArea[5]);
                     int i;for(i=0 ; commArea[i]>=0 ; i+=2)
                         unifyListNoSync(commArea[i], commArea[i+1], 0);
                     //printf("[%d %d]   i/2=%d; commCounter=%d\n", commArea[i-2], commArea[i-1], i/2, commCounter);
                     free(commArea);
                 }
                 printf("Done.\n");
             }
             printf("Ready %d\n", omp_get_thread_num());
         }
         long long num_unify = 0;
         long long load_unify = 0;
         locct[par] = (long long*)callocA64(n*sizeof(long long));
         memory[par] = (int *)mallocA64(n*minPts*sizeof(int));
         double temp[8]__attribute__((aligned(64)));
         printf("thread id %d\n", omp_get_thread_num());
         int commSize = 1024*1024;
         int * commArea = (int *) malloc ((commSize * 2 + 2) * sizeof(int));
         int commCount = 0;
 //        if(par==NUM_THREADS) for(int i=0 ; i<NUM_THREADS ; i++){
 //            omp_set_lock(&testlock);
 //            printf("%2d Testlock granted task_id=%d\n", i, omp_get_thread_num());
 //        }

     EGO_LOOP_TRAN{
         load_unify++;
         BLA2(sum1, 0);
         BLA2(sum2, 1);
         BLA2(sum3, 2);
         BLA2(sum4, 3);
         BLA2(sum5, 4);
         BLA2(sum6, 5);
         BLA2(sum7, 6);
         BLA2(sum8, 7);
     }
     //EGO_CONSOLIDATE{
     }}} FGF_HILBERT_END(i, j); }{
         int curload=0;
         for(int i=loadstart[par] ; i<loadstart[par+1] ; i++)
             for(int s=0 ; s<EGO_stripes ; s++)
                 curload += upper[s][i+nn/8] - lower[s][i+nn/8];
         commArea[2*commCount] = -1 ;
         submitUnifyList(commArea, commField, &commCounter, testlock);
         printf("Consolidate %6.2f %d %d %d %d %d %ld %ld\n",timestamp()-starttimestamp, par, omp_get_thread_num(), loadstart[par], loadstart[par+1]-loadstart[par], curload, num_unify, load_unify);
 #pragma omp critical
         if(++done>=NUM_THREADS)
             omp_unset_lock(&testlock);
         int segments_added=0;
         int curseg = 2*(omp_get_thread_num()%NUM_THREADS) ;
         char * segment=(char*) calloc(2*NUM_THREADS,1);
         while(segments_added < 2 * NUM_THREADS){
             if(!segment[curseg] && omp_test_lock(locks+curseg)){
                 for(int i=curseg*n/NUM_THREADS/16*8 ; i<(curseg+1)*n/NUM_THREADS/16*8 ; i+=8)
                     _mm512_store_si512(counter+i, _mm512_add_epi64(_mm512_load_si512(counter+i),_mm512_load_si512(locct[par]+i)));
                 omp_unset_lock(locks+curseg);
                 segments_added++;
                 segment[curseg] = 1;
             }
             if(++curseg >= 2*NUM_THREADS)
                 curseg = 0;
         }
     }
 }
 }
     }
     printf("This is the point after for par thread=%d!!!!\n", omp_get_thread_num());
 #pragma omp taskwait
     transpose_dx8(EGO_n, EGO_d, EGO_array);}
     printf("Join ready. Cleanup.\n");
     done = NUM_THREADS;
     omp_unset_lock(&testlock);
 //    for(int par=0 ; par<NUM_THREADS ; par++)
 //    #pragma omp task
 //    {
 //        for(int i=0 ; i<n ; i++)
 //            for(int j=0 ; j<minPts && j<locct[par][i] ; j++)
 //                unifyListDbscan(i, memory[par][i*minPts+j], counter[i]>=minPts, counter[memory[par][i*minPts+j]]>=minPts);
 //    }
     for (int i = 0; i < 20; i++)printf("%ld ", counter[i]);
     printf("\n");
     for (int i = n - 20; i < n; i++)printf("%ld ", counter[i]);
     printf("\n");

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
                 if(listArray[j]->minid == i && counter[j] >= minPts) ct++ ;
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
     }
 #ifndef XFOLD
 #define XFOLD 1
 #endif

 void generateSinData(int n, int d, int k, double noise, double* array){
     double rot[k][2][d];
     double trans[k][d];
     for (int i=0 ; i<k*2*d ; i++)
         rot[0][0][i] = drand48() + drand48() + drand48() - 1.5;
     for(int i=0 ; i<k ; i++){
         double h=0;
         for (int j=0 ; j<d ; j++) h+=rot[i][0][j]*rot[i][0][j];
         h = sqrt(h) ;
         for (int j=0 ; j<d ; j++) rot[i][0][j]/=h;
         h=0;
         for (int j=0 ; j<d ; j++) h += rot[i][0][j]*rot[i][1][j];
         for (int j=0 ; j<d ; j++) rot[i][1][j] -= h * rot[i][0][j];
         h=0;
         for (int j=0 ; j<d ; j++) h+=rot[i][1][j]*rot[i][1][j];
         h = sqrt(h) ;
         for (int j=0 ; j<d ; j++) rot[i][1][j]/=h;
     }
     for (int i=0 ; i<k*d ; i++)
         trans[0][i] = (drand48() + drand48() + drand48() - 1.5) * 3.;
     for(int i=0 ; i<n ; i++){
         if(drand48()<noise){
             for (int j=0 ; j<d ; j++)
                 array[i*d+j] = drand48() + drand48() + drand48() - 1.5;
         } else {
             int c=drand48()*k;
             double x=drand48()*2.0-1.0;
             double y=sin(3.14*x);
             for (int j=0 ; j<d ; j++)
                 array[i*d+j] = rot[c][0][j] * x + rot[c][1][j] * y + trans[c][j] + drand48() / 3 / d;
         }
 //        for(int j=0 ; j<d ; j++)
 //            printf("%6.3f ", array[i*d+j]);
 //        printf("\n");
     }
 }
 void generateSinDataGeneral(int n, int d, int k, int cdim, double noise, double* array){
     double rot[k][cdim][d];
     double trans[k][d];
     for (int i=0 ; i<k*cdim*d ; i++)
         rot[0][0][i] = drand48() + drand48() + drand48() - 1.5;
     double sproducts[cdim];
     for(int i=0 ; i<k ; i++){
         for(int a=0 ; a<cdim ; a++){
             for(int b=0 ; b<a ; b++){
                 sproducts[b] = 0;
                 for (int j=0 ; j<d ; j++)
                     sproducts[b] += rot[i][a][j]*rot[i][b][j];
             }
             for(int b=0 ; b<a ; b++)
                 for (int j=0 ; j<d ; j++)
                     rot[i][a][j] -= sproducts[b] * rot[i][b][j];
             double h=0;
             for (int j=0 ; j<d ; j++) h+=rot[i][a][j]*rot[i][a][j];
             h = sqrt(h) ;
             for (int j=0 ; j<d ; j++) rot[i][a][j]/=h;
         }
     }
     for (int i=0 ; i<k*d ; i++)
         trans[0][i] = (drand48() + drand48() + drand48() - 1.5) * 3.;
     for(int i=0 ; i<n ; i++){
         if(drand48()<noise){
             for (int j=0 ; j<d ; j++)
                 array[i*d+j] = drand48() + drand48() + drand48() - 1.5;
         } else {
             int c=drand48()*k;
             sproducts[0] = 1;
             for(int j=1 ; j<cdim ; j++){
                 sproducts[j] = drand48()*2.0-1.0;
                 sproducts[0] *= sin(3.14*sproducts[j]);
             }
             for (int j=0 ; j<d ; j++){
                 array[i*d+j] = trans[c][j] + drand48() / 10 / d;
                 for(int a=0 ; a<cdim ; a++)
                     array[i*d+j] += rot[c][a][j] * sproducts[a];
             }
         }
         for(int j=0 ; j<d ; j++)
             printf("%6.3f ", array[i*d+j]);
         printf("\n");
     }
 }
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
                 int allelements = seedcounter; // elements in seedlist plus non-waiting threads
                 cores = 1;
 #pragma omp parallel for
                 for (int par = 0; par < NUM_THREADS; par++) {
                     int ii;
                     while (allelements > 0) {
                         // WAIT for the counting semaphor
                         omp_set_lock(&delay);
 #pragma omp critical(mutex)
                         {
                             if (allelements > 0) {
                                 ii = seedlist[--seedcounter];
                                 if (seedcounter > 0)
                                     omp_unset_lock(&delay);
                             }
                         }
                         if(allelements <= 0) break;
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
                                                     allelements ++;
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
                                         allelements ++;
                                         if (seedcounter == 1) omp_unset_lock(&delay);
                                     }
                             }
                             allelements--;
                         }
                     }
                     omp_unset_lock(&delay);
                 }
                 if (cores > 500)
                     printf("Cluster %d p=%d crs=%d\n", curId, p, cores);
                 curId++;
             }
         }
     }
 }

 void warshallSimple(int d, int *A){
     for (int k=0 ; k<d ; k++)
         for(int i=0; i<d ; i++) {
             if(A[d*i+k])
                 for (int j=0 ; j<d ; j++)
 //                    A[d*i+j]|=A[d*k+j];
                     if(A[d*k+j])
                         A[d*i+j] = 1;
         }
 }

 void warshallSimpleOMP(int d, int *A){
     for (int k=0 ; k<d ; k++)
 # pragma omp parallel for
         for (int par = 0; par < NUM_THREADS; par++) {
             for (int i = par * d / NUM_THREADS ; i < (par + 1) * d / NUM_THREADS ; i++)
             if(A[d*i+k])
                 for (int j=0 ; j<d ; j++)
 //                    A[d*i+j]|=A[d*k+j];
                     if(A[d*k+j])
                         A[d*i+j] = 1;
         }
 }

 void pack8(int d, int * A, unsigned char * B){
     int dd=(d+7)/8;
     for(int i=0 ; i<d ; i++){
         for(int j=0; j<dd ; j++){
             int * AA=A+i*d+j*8;
 //            B[i*dd+j]=AA[0]*128+AA[1]*64+AA[2]*32+AA[3]*16+AA[4]*8+AA[5]*4+AA[6]*2+AA[7];
             B[i*dd+j]=((AA[0]&1)<<7) | ((AA[1]&1)<<6) | ((AA[2]&1)<<5) | ((AA[3]&1)<<4) | ((AA[4]&1)<<3) | ((AA[5]&1)<<2) | ((AA[6]&1)<<1) | ((AA[7]&1));
         }
         B[i*dd+dd-1] &= 255<<((-d)&7) ;
     }
 }
 void warshallSimpleAVX(int d, unsigned char *A){
     int dd = (d+7)/8 ;
     unsigned char mask [32] ;
     for(int i=0 ; i<d/8%32 ; i++) mask[i] = 255;
     for(int i=d/8%32 ; i<32 ; i++) mask[i] = 0;
     mask[d/8%32] = 255<<((-d)&7) ;
     __m256i maskr = _mm256_loadu_si256((__m256i*)mask);
 //        for(int j=0; j<32 ; j++)
 //            printf("%02x ", mask[j]);
 //        printf("\n");
     for (int k=0 ; k<d ; k++){
         for (int i=0 ; i<d ; i++){
 //            printf("%d",(A[dd*i+k/8]>>((7-k)&7))&1);
             if((A[dd*i+k/8]>>((7-k)&7))&1){
                 for (int j=0 ; j<dd ; j+=32){
                     __m256i h = _mm256_loadu_si256 ((__m256i*)(A + dd*k + j));
                     if(j+32>=dd)
                         h = _mm256_and_si256(h, maskr);
                     h = _mm256_or_si256 (h, _mm256_loadu_si256 ((__m256i*)(A + dd*i + j)));
                     _mm256_storeu_si256 ((__m256i*)(A + dd*i + j), h);
                 }
             }
         }
 //        printf("\n");
    }
 }
 void warshallSimpleAVXOMP(int d, unsigned char *A){
     int dd = (d+7)/8 ;
     unsigned char mask [32] ;
     for(int i=0 ; i<d/8%32 ; i++) mask[i] = 255;
     for(int i=d/8%32 ; i<32 ; i++) mask[i] = 0;
     mask[d/8%32] = 255<<((-d)&7) ;
 //        for(int j=0; j<32 ; j++)
 //            printf("%02x ", mask[j]);
 //        printf("\n");
     for (int k = 0; k < d; k++) {
 #pragma omp parallel for
         for (int par = 0; par < NUM_THREADS; par++) {
             __m256i maskr = _mm256_loadu_si256((__m256i*) mask);
             for (int i = (par + (k & 1)) * d / NUM_THREADS - (k & 1); i < (par + 1) * d / NUM_THREADS && i >= par*d/NUM_THREADS ; i -= (k & 1)*2 - 1) {
                 // printf("%d",(A[dd*i+k/8]>>((7-k)&7))&1);
                 // printf("%d %d\n", k, i);
                 if ((A[dd * i + k / 8]>>((7 - k)&7))&1) {
                     for (int j = 0; j < dd; j += 32) {
                         __m256i h = _mm256_loadu_si256((__m256i*) (A + dd * k + j));
                         if (j + 32 >= dd)
                             h = _mm256_and_si256(h, maskr);
                         h = _mm256_or_si256(h, _mm256_loadu_si256((__m256i*) (A + dd * i + j)));
                         _mm256_storeu_si256((__m256i*) (A + dd * i + j), h);
                     }
                 }
             }
         }
         //printf("\n");
     }
 }

 void warshallSimpleScramble(int d, unsigned char *A){
     int primzahlen[] = {
    61,    67,    71,    73,    79,    83,    89,    97,   101,   103,   107,   109,   113,   127,   131,   137,
   139,   149,   151,   157,   163,   167,   173,   179,   181,   191,   193,   197,   199,   211,   223,   227,
   229,   233,   239,   241,   251,   257,   263,   269,   271,   277,   281,   283,   293,   307,   311,   313,
   317,   331,   337,   347,   349,   353,   359,   367,   373,   379,   383,   389,   397,   401,   409,   419};
     int dd = (d+7)/8 ;
     unsigned char mask [32] ;
     for(int i=0 ; i<d/8%32 ; i++) mask[i] = 255;
     for(int i=d/8%32 ; i<32 ; i++) mask[i] = 0;
     mask[d/8%32] = 255<<((-d)&7) ;
     for (int k = 0; k < d; k++) {
 #pragma omp parallel for
         for (int pp = 0; pp < NUM_THREADS; pp++) {
             int par = pp*primzahlen[k%64] % NUM_THREADS ;
             __m256i maskr = _mm256_loadu_si256((__m256i*) mask);
             for (int i = (par + (k & 1)) * d / NUM_THREADS - (k & 1); i < (par + 1) * d / NUM_THREADS && i >= par*d/NUM_THREADS ; i -= (k & 1)*2 - 1) {
                 // printf("%d",(A[dd*i+k/8]>>((7-k)&7))&1);
                 // printf("%d %d\n", k, i);
                 if ((A[dd * i + k / 8]>>((7 - k)&7))&1) {
                     for (int j = 0; j < dd; j += 32) {
                         __m256i h = _mm256_loadu_si256((__m256i*) (A + dd * k + j));
                         if (j + 32 >= dd)
                             h = _mm256_and_si256(h, maskr);
                         h = _mm256_or_si256(h, _mm256_loadu_si256((__m256i*) (A + dd * i + j)));
                         _mm256_storeu_si256((__m256i*) (A + dd * i + j), h);
                     }
                 }
             }
         }
         //printf("\n");
     }
 }

 void warshallSimpleBarrier(int d, unsigned char *A) {
     int dd = (d + 7) / 8;
     unsigned char mask [32];
     for (int i = 0; i < d / 8 % 32; i++) mask[i] = 255;
     for (int i = d / 8 % 32; i < 32; i++) mask[i] = 0;
     mask[d / 8 % 32] = 255 << ((-d)&7);
     int par;
 #pragma omp parallel private(par)
     {
 #pragma omp for
     for (int par = 0; par < NUM_THREADS; par++) {
         __m256i maskr = _mm256_loadu_si256((__m256i*) mask);
         for (int k = 0; k < d; k++) {
 //#pragma omp barrier
             for (int i = (par + (k & 1)) * d / NUM_THREADS - (k & 1); i < (par + 1) * d / NUM_THREADS && i >= par * d / NUM_THREADS; i -= (k & 1)*2 - 1) {
                 if ((A[dd * i + k / 8]>>((7 - k)&7))&1) {
                     for (int j = 0; j < dd; j += 32) {
                         __m256i h = _mm256_loadu_si256((__m256i*) (A + dd * k + j));
                         if (j + 32 >= dd)
                             h = _mm256_and_si256(h, maskr);
                         h = _mm256_or_si256(h, _mm256_loadu_si256((__m256i*) (A + dd * i + j)));
                         _mm256_storeu_si256((__m256i*) (A + dd * i + j), h);
                     }
                 }
             }
         }
 //#pragma omp barrier
     }}
 }

 void choleskyNonrekBlockoutput(int d) {
     for (int a=0 ; a<d ; a+=1) {
         if (a) {
             int b=(a & -a) ; // size of Block: greatest power of 2 that divides a with no rest
             // Block (a .. min(a+b,d), a-b .. a)
             int c = a+b>d ? d : a+b ;
             printf("%.1f, %.1f\n%.1f, %.1f\n%.1f, %.1f\n%.1f, %.1f\n%.1f, %.1f\n\n", a-b-0.5, a-0.5, a-b-0.5, c-0.5, a-0.5, c-0.5, a-0.5, a-0.5, a-b-0.5, a-0.5) ;
         }
         // Triangle (a..a+4, a..a+4)
 //        printf("%d, %d\n%d, %d\n%d, %d\n%d, %d\n\n", a, a, a+4, a, a+4, a+4, a, a) ;
     }
 }

 void warshallHilloop(int d, int *A) {
     for (int a=0 ; a<d ; a+=4) {
         if (a) {
             int b=(a & -a) ; // size of Block: greatest power of 2 that divides a with no rest
             // Block (a .. min(a+b,d), a-b .. a)
             int c = a+b>d ? d : a+b ;
             int kk, i ;
             FUR_HILBERT_START(kk, i, 0, b, a, c) {
                 int k = a-kk-1 ;
                   printf("%d %d\n", k, i);
               if(A[d*i+k])
                     for (int j=0 ; j<d ; j++)
                         A[d*i+j] |= A[d*k+j] ;
             } FUR_HILBERT_END(kk, i) ;
         }
         int i=a+4 ;
         for (int k=a ; k<a+4 ; k++)
             for(i+=(k&1)*2-1 ; i>=0 && i<a+4 ; i+=(k&1)*2-1){
                   printf("%d %d\n", k, i);
               if(A[d*i+k])
                     for (int j=0 ; j<d ; j++)
                         A[d*i+j] |= A[d*k+j] ;
             }
     }
 }

 void warshallHilloopTwoWay(int d, int *A, int v) {
 //    int v = 0;
     int steps = 2;
     for (int k = 0; k < steps; k++)
         for (int i = 0; i < steps; i++)
             if (i != k) {
                 if (v) printf("%d %d\n", k, i);
                 if (A[d * i + k]) for (int j = 0; j < d; j++) A[d * i + j] |= A[d * k + j];
             }
     if(v) printf("%d %d\n",0,0);

     for (int a = steps; a < d; a += steps) {
         int b = (a & -a); // size of Block: greatest power of 2 that divides a with no rest
         // Block (a .. min(a+b,d), a-b .. a)
         int c = a + b > d ? d : a + b;
         int kk, i;

         FUR_HILBERT_START(kk, i, 0, b, a, c) {
             int k = a - kk - 1;
             if (v) printf("%d %d\n", k, i);
             if (A[d * i + k])
                 for (int j = 0; j < d; j++)
                     A[d * i + j] |= A[d * k + j];
         }
         FUR_HILBERT_END(kk, i);


 //        CANO2LOOP_START(i, kk, a - b, a, a, c) {
 //            int k = kk;
 //            if (v) printf("%d %d\n", k, i);
 //            if (A[d * i + k])
 //                for (int j = 0; j < d; j++)
 //                    A[d * i + j] |= A[d * k + j];
 //        }
 //        CANO2LOOP_END(i, kk);

 //        if(v) printf("%d %d\n",a,a);
         for (int k = a; k < a + steps && k < d ; k++)
             for (int i = (a+steps-1)*(k&1); i < a + steps && i < d && i>=0 ; i-=(k&1)*2-1)
                 if (i != k) {
                     if (v) printf("%d %d\n", k, i);
                     if (A[d * i + k]) for (int j = 0; j < d; j++) A[d * i + j] |= A[d * k + j];
                 }
         //        k=a; i=a+1; if(v) printf("%d %d\n", k, i); if (A[d * i + k]) for (int j = 0; j < d; j++) A[d * i + j] |= A[d * k + j];
         //        k=a+1; i=a; if(v) printf("%d %d\n", k, i); if (A[d * i + k]) for (int j = 0; j < d; j++) A[d * i + j] |= A[d * k + j];
         //              if(A[d*i+k])
         //                    for (int j=0 ; j<d ; j++)
         //                        A[d*i+j] |= A[d*k+j] ;
     }
 }

 void warshallHilloopAVXOMP(int d, unsigned char *A) {
     int dd = (d + 7) / 8;
     unsigned char mask [32];
     for (int i = 0; i < d / 8 % 32; i++) mask[i] = 255;
     for (int i = d / 8 % 32; i < 32; i++) mask[i] = 0;
     mask[d / 8 % 32] = 255 << ((-d)&7);
     for (int a = 0; a < d; a += 4) {
         if (a) {
             int b = (a & -a); // size of Block: greatest power of 2 that divides a with no rest
             // Block (a .. min(a+b,d), a-b .. a)
             int c = a + b > d ? d : a + b;
             int ca = c - a;
             int num_threads = min(NUM_THREADS, ca);
 #pragma omp parallel for
             for (int par = 0; par < num_threads; par++) {
                 __m256i maskr = _mm256_loadu_si256((__m256i*) mask);
                 int kk, i;

                 FUR_HILBERT_START(kk, i, 0, b, a + par * ca / num_threads, a + (par + 1) * ca / num_threads) {
                     int k = a - kk - 1;
 //                    printf("%d %d %d\n", par, k, i);
                     if ((A[dd * i + k / 8]>>((7 - k)&7))&1) {
                         for (int j = 0; j < dd; j += 32) {
                             __m256i h = _mm256_loadu_si256((__m256i*) (A + dd * k + j));
                             if (j + 32 >= dd)
                                 h = _mm256_and_si256(h, maskr);
                             h = _mm256_or_si256(h, _mm256_loadu_si256((__m256i*) (A + dd * i + j)));
                             _mm256_storeu_si256((__m256i*) (A + dd * i + j), h);
                         }
                     }
                 }
                 FUR_HILBERT_END(kk, i);
             }
         }
         for (int k = a; k < a + 4; k++) {
             int a4 = a + 4;
             int num_threads = min(NUM_THREADS, a4);
 #pragma omp parallel for
             for (int par = 0; par < num_threads; par++) {
                 __m256i maskr = _mm256_loadu_si256((__m256i*) mask);
                 for (int i = (par + (~k & 1)) * a4 /num_threads - (~k & 1) ; i >= par * a4 / num_threads && i < (par + 1) * a4 / num_threads; i += (k & 1)*2 - 1) {
 //                    printf("%d %d %d\n", par, k, i);
                     if ((A[dd * i + k / 8]>>((7 - k)&7))&1) {
                         for (int j = 0; j < dd; j += 32) {
                             __m256i h = _mm256_loadu_si256((__m256i*) (A + dd * k + j));
                             if (j + 32 >= dd)
                                 h = _mm256_and_si256(h, maskr);
                             h = _mm256_or_si256(h, _mm256_loadu_si256((__m256i*) (A + dd * i + j)));
                             _mm256_storeu_si256((__m256i*) (A + dd * i + j), h);
                         }
                     }
                 }
             }
         }
     }
 }

 void warshallHilloopTest(int d, int *A) {
     for (int a = 0; a < d; a += 4) {
         int i = a + 4;
         for (int k = a; k < a + 4; k++)
             for (i += (k & 1)*2 - 1; i >= 0 && i < a + 4; i += (k & 1)*2 - 1) {
                 //                printf("%d %d\n", k, i);
                 if (A[d * i + k])
                     for (int j = 0; j < d; j++)
                         A[d * i + j] |= A[d * k + j];
             }
     }
     for (int a = 0; a < d; a += 4) {
         if (a) {
             int b = (a & -a); // size of Block: greatest power of 2 that divides a with no rest
             // Block (a .. min(a+b,d), a-b .. a)
             int c = a + b > d ? d : a + b;
             int kk, i;

             FUR_HILBERT_START(kk, i, 0, b, a, c) {
                 int k = a - kk - 1;
                 //                printf("%d %d\n", k, i);
                 if (A[d * i + k])
                     for (int j = 0; j < d; j++)
                         A[d * i + j] |= A[d * k + j];
             }
             FUR_HILBERT_END(kk, i);
         }
     }
 }
 void testWarshall(){
     int n = 10000;
     int n8 = (n+7)/8;
     int *testA = (int *)calloc(n*n,sizeof(int));
     int *testB = (int *)calloc(n*n,sizeof(int));
     unsigned char *testC = (unsigned char *) calloc(n*n8, sizeof(char));
     unsigned char *testD = (unsigned char *) calloc(n*n8, sizeof(char));
     for(int i=1; i<n ; i++)
         for(int j=0; j<i ; j++)
             if(i%3 == j%3 && drand48()<0.0005)
                 testA[i*n+j] = testA[j*n+i] = testB[i*n+j] = testB[j*n+i] = 1;
     for(int i=0; i<100 ; i++){
         for(int j=0; j<100 ; j++)
             if(testA[i*n+j]) printf("1");
             else printf(" ");
         printf("\n");
     }
     warshallHilloop(23, testA);
     //choleskyNonrekBlockoutput(23);
     exit(0);
     pack8(n, testA, testC);
     pack8(n, testA, testD);
     for(int i=0; i<100 ; i++){
         for(int j=0; j<13 ; j++)
             printf("%02x ", testC[i*n8+j]);
         printf("\n");
     }
     stop(); stopc();
     warshallSimpleAVXOMP(n,testD);
     printf("warshallSimpleAVXOMP   | %8.3f | %8.3f\n", stopc(), stop());
 //    warshallSimpleScramble(10000, testC);
 //    printf("warshallSimpleScramble  | %8.3f | %8.3f\n", stopc(), stop());
     warshallHilloopAVXOMP(n, testC);
     printf("warshallHilloopAVXOMP  | %8.3f | %8.3f\n", stopc(), stop());
     warshallSimpleOMP(n,testB);
     printf("warshallSimpleOMP      | %8.3f | %8.3f\n", stopc(), stop());
     warshallHilloopTwoWay(n,testA, 0);
     printf("warshallHilloopTwoWay  | %8.3f | %8.3f\n", stopc(), stop());
     for(int i=0; i<100 ; i++){
         for(int j=0; j<100 ; j++)
             if(testA[i*n+j]==testB[i*n+j])
                 if(testA[i*n+j]) printf("|");
                 else printf(" ");
             else printf("X");
         printf("\n");
     }
     for(int i=0; i<100 ; i++){
         for(int j=0; j<13 ; j++){
             if(testC[i*n8+j]==testD[i*n8+j]) printf(" "); else printf("X %02x ", testD[i*n8+j]);
             printf("%02x ", testC[i*n8+j]);
         }
         printf("\n");
     }
     exit(0);
 }

 double entropy(double *d, int n, double s){
     double denom = 0 ;
 //    s *= 2 ; // THIS IS WRONG !!!
     for(int i=0 ; i<n ; i++)
         denom += exp(-d[i]/s);
     double H = 0;
     for(int i=0 ; i<n ; i++){
         double pij = exp(-d[i]/s)/denom ;
 //        printf("%f ",pij);
         if(pij>0)
             H -= pij * log(pij) ;
     }
     return H/log(2.);
 }
 double perpsearch(double *d, int n, double perplexity) {
     // d[0] ... d[n-1]: squared euclidean distances of nearest neighbors
     double l2 = log(2.) ;
     perplexity = log(perplexity)/l2 ;
     printf("goal: %f\n",perplexity);
     double s = 0 ;
 //    for(int i=0 ; i<n ; i++)
 //        if(s == 0 || d[i] < s)
 //            s = d[i] ;
     double lsl = -5 ;
     double lpl = entropy(d,n,exp(lsl));
 //    for(int i=0 ; i<n ; i++)
 //        if(d[i] > s)
 //            s = d[i] ;
 //    if(s<=0) return 0 ;
     double lsh = 5 ;
     double lph = entropy(d,n,exp(lsh));
     for(int i=0 ;; i++){
         s = (lsh-lsl)/(lph-lpl)*(perplexity-lpl)+lsl;
         // s = (lsl+lsh)/2.; //0.999*s +
         double lp = entropy(d,n,exp(s));
         printf("%d %f %f %f %f %f %f\n", i, lsl, s, lsh, lpl, lp, lph);
         if(abs(lp-perplexity)/perplexity < 1e-10 || i>50)
             return exp(s) ;
         if(lp > perplexity){
             lph = lp ;
             lsh = s ;
         } else {
             lpl = lp ;
             lsl = s ;
         }
     }
 }
 void poweriter(int n, int d, double *a, double *eigenv){
     int num_iter = 10 ;
     double mean[d] ;
     double cov[d][d] ;
     for (int j=0 ; j<d ; j++){
         mean[j] = a[j] ;
         for (int i=1 ; i<n ; i++) mean[j] += a[d*i+j];
         mean[j] /= n ;
     }
     for (int i=0 ; i<d ; i++)
         for (int j=0 ; j<d ; j++){
             cov[i][j] = (a[i]-mean[i])*(a[j]-mean[j]) ;
             for (int k=1 ; k<n ; k++)
                 cov[i][j] += (a[k*d+i]-mean[i])*(a[k*d+j]-mean[j]) ;
             cov[i][j] /= n ;
             printf(j==d-1?"%f\n":"%f ", cov[i][j]);
         }
     for (int i=0 ; i<d ; i++)
         eigenv[i] = drand48() ;
     double h[d] ;
     for (int i=0 ; i<num_iter ; i++){
         double length = 0.;
         for(int j=0 ; j<d ; j++) {
             h[j] = eigenv[0]*cov[j][0] ;
             for (int k=1 ; k<d ; k++)
                 h[j] += eigenv[k]*cov[j][k] ;
             length += h[j]*h[j] ;
         }
         length=sqrt(length);
         for(int j=0 ; j<d ; j++)
             eigenv[j] = h[j] / length ;
         printf("%f %f %f %f\n", eigenv[0], eigenv[1], eigenv[2], eigenv[3]);
     }
 }
 double iris[] = {
     5.1, 3.5, 1.4, 0.2, 4.9, 3.0, 1.4, 0.2, 4.7, 3.2, 1.3, 0.2, 4.6, 3.1, 1.5, 0.2,
     5.0, 3.6, 1.4, 0.2, 5.4, 3.9, 1.7, 0.4, 4.6, 3.4, 1.4, 0.3, 5.0, 3.4, 1.5, 0.2,
     4.4, 2.9, 1.4, 0.2, 4.9, 3.1, 1.5, 0.1, 5.4, 3.7, 1.5, 0.2, 4.8, 3.4, 1.6, 0.2,
     4.8, 3.0, 1.4, 0.1, 4.3, 3.0, 1.1, 0.1, 5.8, 4.0, 1.2, 0.2, 5.7, 4.4, 1.5, 0.4,
     5.4, 3.9, 1.3, 0.4, 5.1, 3.5, 1.4, 0.3, 5.7, 3.8, 1.7, 0.3, 5.1, 3.8, 1.5, 0.3,
     5.4, 3.4, 1.7, 0.2, 5.1, 3.7, 1.5, 0.4, 4.6, 3.6, 1.0, 0.2, 5.1, 3.3, 1.7, 0.5,
     4.8, 3.4, 1.9, 0.2, 5.0, 3.0, 1.6, 0.2, 5.0, 3.4, 1.6, 0.4, 5.2, 3.5, 1.5, 0.2,
     5.2, 3.4, 1.4, 0.2, 4.7, 3.2, 1.6, 0.2, 4.8, 3.1, 1.6, 0.2, 5.4, 3.4, 1.5, 0.4,
     5.2, 4.1, 1.5, 0.1, 5.5, 4.2, 1.4, 0.2, 4.9, 3.1, 1.5, 0.2, 5.0, 3.2, 1.2, 0.2,
     5.5, 3.5, 1.3, 0.2, 4.9, 3.6, 1.4, 0.1, 4.4, 3.0, 1.3, 0.2, 5.1, 3.4, 1.5, 0.2,
     5.0, 3.5, 1.3, 0.3, 4.5, 2.3, 1.3, 0.3, 4.4, 3.2, 1.3, 0.2, 5.0, 3.5, 1.6, 0.6,
     5.1, 3.8, 1.9, 0.4, 4.8, 3.0, 1.4, 0.3, 5.1, 3.8, 1.6, 0.2, 4.6, 3.2, 1.4, 0.2,
     5.3, 3.7, 1.5, 0.2, 5.0, 3.3, 1.4, 0.2, 7.0, 3.2, 4.7, 1.4, 6.4, 3.2, 4.5, 1.5,
     6.9, 3.1, 4.9, 1.5, 5.5, 2.3, 4.0, 1.3, 6.5, 2.8, 4.6, 1.5, 5.7, 2.8, 4.5, 1.3,
     6.3, 3.3, 4.7, 1.6, 4.9, 2.4, 3.3, 1.0, 6.6, 2.9, 4.6, 1.3, 5.2, 2.7, 3.9, 1.4,
     5.0, 2.0, 3.5, 1.0, 5.9, 3.0, 4.2, 1.5, 6.0, 2.2, 4.0, 1.0, 6.1, 2.9, 4.7, 1.4,
     5.6, 2.9, 3.6, 1.3, 6.7, 3.1, 4.4, 1.4, 5.6, 3.0, 4.5, 1.5, 5.8, 2.7, 4.1, 1.0,
     6.2, 2.2, 4.5, 1.5, 5.6, 2.5, 3.9, 1.1, 5.9, 3.2, 4.8, 1.8, 6.1, 2.8, 4.0, 1.3,
     6.3, 2.5, 4.9, 1.5, 6.1, 2.8, 4.7, 1.2, 6.4, 2.9, 4.3, 1.3, 6.6, 3.0, 4.4, 1.4,
     6.8, 2.8, 4.8, 1.4, 6.7, 3.0, 5.0, 1.7, 6.0, 2.9, 4.5, 1.5, 5.7, 2.6, 3.5, 1.0,
     5.5, 2.4, 3.8, 1.1, 5.5, 2.4, 3.7, 1.0, 5.8, 2.7, 3.9, 1.2, 6.0, 2.7, 5.1, 1.6,
     5.4, 3.0, 4.5, 1.5, 6.0, 3.4, 4.5, 1.6, 6.7, 3.1, 4.7, 1.5, 6.3, 2.3, 4.4, 1.3,
     5.6, 3.0, 4.1, 1.3, 5.5, 2.5, 4.0, 1.3, 5.5, 2.6, 4.4, 1.2, 6.1, 3.0, 4.6, 1.4,
     5.8, 2.6, 4.0, 1.2, 5.0, 2.3, 3.3, 1.0, 5.6, 2.7, 4.2, 1.3, 5.7, 3.0, 4.2, 1.2,
     5.7, 2.9, 4.2, 1.3, 6.2, 2.9, 4.3, 1.3, 5.1, 2.5, 3.0, 1.1, 5.7, 2.8, 4.1, 1.3,
     6.3, 3.3, 6.0, 2.5, 5.8, 2.7, 5.1, 1.9, 7.1, 3.0, 5.9, 2.1, 6.3, 2.9, 5.6, 1.8,
     6.5, 3.0, 5.8, 2.2, 7.6, 3.0, 6.6, 2.1, 4.9, 2.5, 4.5, 1.7, 7.3, 2.9, 6.3, 1.8,
     6.7, 2.5, 5.8, 1.8, 7.2, 3.6, 6.1, 2.5, 6.5, 3.2, 5.1, 2.0, 6.4, 2.7, 5.3, 1.9,
     6.8, 3.0, 5.5, 2.1, 5.7, 2.5, 5.0, 2.0, 5.8, 2.8, 5.1, 2.4, 6.4, 3.2, 5.3, 2.3,
     6.5, 3.0, 5.5, 1.8, 7.7, 3.8, 6.7, 2.2, 7.7, 2.6, 6.9, 2.3, 6.0, 2.2, 5.0, 1.5,
     6.9, 3.2, 5.7, 2.3, 5.6, 2.8, 4.9, 2.0, 7.7, 2.8, 6.7, 2.0, 6.3, 2.7, 4.9, 1.8,
     6.7, 3.3, 5.7, 2.1, 7.2, 3.2, 6.0, 1.8, 6.2, 2.8, 4.8, 1.8, 6.1, 3.0, 4.9, 1.8,
     6.4, 2.8, 5.6, 2.1, 7.2, 3.0, 5.8, 1.6, 7.4, 2.8, 6.1, 1.9, 7.9, 3.8, 6.4, 2.0,
     6.4, 2.8, 5.6, 2.2, 6.3, 2.8, 5.1, 1.5, 6.1, 2.6, 5.6, 1.4, 7.7, 3.0, 6.1, 2.3,
     6.3, 3.4, 5.6, 2.4, 6.4, 3.1, 5.5, 1.8, 6.0, 3.0, 4.8, 1.8, 6.9, 3.1, 5.4, 2.1,
     6.7, 3.1, 5.6, 2.4, 6.9, 3.1, 5.1, 2.3, 5.8, 2.7, 5.1, 1.9, 6.8, 3.2, 5.9, 2.3,
     6.7, 3.3, 5.7, 2.5, 6.7, 3.0, 5.2, 2.3, 6.3, 2.5, 5.0, 1.9, 6.5, 3.0, 5.2, 2.0,
     6.2, 3.4, 5.4, 2.3, 5.9, 3.0, 5.1, 1.8
 };

 double irissig [150];
 double eigenv [4];

 void fillperplex (int n, int d, double *a, double *s, double perplexity){
     double dists[n];
     for(int i=0 ; i<n ; i++){
         for(int j=0 ; j<n ; j++) {
             dists[j-(j>i)] = (a[i*d]-a[j*d])*(a[i*d]-a[j*d]) ;
             for(int k=1 ; k<d ; k++)
                 dists[j-(j>i)] += (a[i*d+k]-a[j*d+k])*(a[i*d+k]-a[j*d+k]) ;
         }
         // TEST TEST TEST
 //        for(int i=0 ; i<149 ; i++) printf("%f ", dists[i]);
 //        printf("\n");
 //        entropy(dists, 149, 1.);
 //        exit(0);
         s[i] = perpsearch(dists, n-1, perplexity);
         printf("\n%d %f %f\n", i, s[i], sqrt(s[i]));
     }
 }
 int PERPJ_d ;
 double * PERPJ_ev ;
 int eigenVectorCompare(const void *a, const void *b) {
     double * A = (double *) a;
     double * B = (double *) b;
     double h = 0.;
     for(int i=0 ; i<PERPJ_d ; i++)
         h += (A[i]-B[i]) * PERPJ_ev[i] ;
     if(h<0) return -1;
     if(h>0) return 1;
     return 0;
 }
 void eigenVectorOrdering(int n, int d, double *ev, double *array) {
     PERPJ_ev = ev;
     PERPJ_d = d;
     min_size_qsort = n/NUM_THREADS*2 ;
 #pragma omp parallel
 #pragma omp master
     omp_qsort(array, n, d * sizeof (double), eigenVectorCompare);
 #pragma omp taskwait
 }
 double * createProjection(int n, int d, double *ev, double *array){
     double *result = (double *) calloc(n, sizeof(double));
     for(int i=0 ; i<n ; i++)
         for (int j=0 ; j<d ; j++)
             result[i]+=array[i*d+j]*ev[j];
     return result ;
 }
 int * perpjFillList(int n, int d, double *proj, double *sig, double factor) {
     int nn = ceilpowtwo(n);
     int *result = (int *) callocA64(nn*2*sizeof(int)) + nn;
 #pragma omp parallel for reduction(max:result[:nn])
     for(int par=0 ; par<NUM_THREADS ; par++){
         int j;
         for(int i=j=par*n/NUM_THREADS ; i<(par+1)*n/NUM_THREADS ; i++){
             while (j<n && proj[j] <= proj[i]+factor*sig[i]) j++;
             while (j>i && proj[j-1] > proj[i]+factor*sig[i]) j--;
             result[i] = j;
             for(int k=i; k>=0 && proj[k] >= proj[i] - factor*sig[i] ; k--)
                 if(i > result[k])
                     result[k] = i;
         }
     }
     result -= nn;
     for (int i = nn - 1; i > 0; i--)
         result[i] = max(result[2 * i], result[2 * i + 1]);
     return result ;
 }
 # define PJ_PARALLEL_TRAN(n,d,array,factor) {\
     int PJ_n = (n);\
     int PJ_d = (d);\
     double PJ_sig = (sig);\
     double * PJ_array = (array);\
     double PJ_factor = (factor);\
     int PJ_blocks = (PJ_d + KBLOCK - 1) / KBLOCK;\
     unsigned long long usedload[NUM_THREADS]; for(int i=0; i<NUM_THREADS; i++) usedload[i]=0ull;\
     omp_lock_t criticallock; omp_init_lock(&criticallock); int scritical = -1;\
     double PJ_ev[PJ_d];\
     poweriter(PJ_n, PJ_d, PJ_array, PJ_ev);\
     eigenVectorOrdering(PJ_n, PJ_d, PJ_ev, PJ_array);\
     int nn = ceilpowtwo(PJ_n);\
     double *self = callocA64(sizeof (double) * PJ_n * PJ_blocks);\
     int * PJ_stripe = perpjFillList(PJ_n, PJ_d, PJ_proj, PJ_sig, PJ_factor);\
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
     unsigned long long *sloadcumcum = (unsigned long long *) malloc(EGO_stripes * (n+7)/8 * sizeof(unsigned long long));\
     sloadcumcum[EGO_stripes-1]=upper[EGO_stripes-1][nn/8]-lower[EGO_stripes-1][nn/8];\
     for(int s=EGO_stripes-1; s>=0 ; s--)\
         sloadcumcum[s] = sloadcumcum[s+1] + upper[s][nn/8] - lower[s][nn/8];\
     for(int i=1; i<(n+7)/8 ; i++)\
         for(int s=EGO_stripes-1; s>=0 ; s--)\
             sloadcumcum[i*EGO_stripes+s] = sloadcumcum[(i-1)*EGO_stripes+s] + upper[s][nn/8+i] - lower[s][nn/8+i];\
     for(int i=0; i<(n+7)/8 ; i++)\
         for(int s=EGO_stripes-2; s>=0 ; s--)\
             sloadcumcum[i*EGO_stripes+s] += sloadcumcum[i*EGO_stripes+s+1];\
     unsigned long long * assignedload = (unsigned long long *) calloc(NUM_THREADS, sizeof (unsigned long long));\
     long long overall_load = 0;\
     for (int i = 0; i < EGO_n / 8; i++)\
         for (int j = 0; j < EGO_stripes; j++)\
             overall_load += upper[j][i + nn / 8] - lower[j][i + nn / 8];\
     int *loadstart = (int *) calloc((NUM_THREADS + 1) * EGO_stripes, sizeof(int));\
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

 long long * costref;
 int cmp_reorder_dim(const void * a, const void *b){
     long long A = costref[*(int*) a];
     long long B = costref[*(int*) b];
     if(A<B) return -1;
     if(A>B) return 1;
     return 0;
 }

 void outputStatistics(int n, int d, double epsilon, double *array) {
     printf("Statistics\n");
     //    int first_int[d], last_int[d];
     //    for (int j=0 ; j<d ; j++)
     //        first_int[j] = last_int[j] = (int) floor(array[j]/epsilon);
     //    for (int i=0 ; i<n ; i++)
     //        for (int j=0 ; j<d ; j++){
     //            first_int[j] = min(first_int[j], (int) floor(array[i*d+j]/epsilon));
     //            last_int[j] = max(last_int[j], (int) floor(array[i*d+j]/epsilon));
     //        }
     costref = (long long *) malloc(d*sizeof(long long));
     for (int j = 0; j < d; j++) {
         double first_int = array[j];
         double last_int = first_int;
         int testcount = 0;
 #pragma omp parallel for reduction(min:first_int) reduction(max:last_int) reduction(+:testcount)
         for(int par=0 ; par<NUM_THREADS ; par++)
         for (int i = par+1; i < n; i+=(int)(-log(1.-drand48())*100.*NUM_THREADS)){
             first_int = min(first_int, array[i * d + j]);
             last_int = max(last_int, array[i * d + j]);
             testcount ++;
         }
         first_int = floor(first_int/epsilon);
         last_int = floor(last_int/epsilon);
         if (first_int < -100000 || last_int > 100000) {
             printf("%3d %e %e ------------------ (VALUE RANGE EXCEEDED)\n", j, first_int, last_int);
         } else {
             int size = (int) last_int - (int) first_int + 1;
             long long * histo = (long long *) calloc(size, sizeof (long long));
             //memset(histo, 0, size*sizeof(long long));
             int testcount2 = 0;
 #pragma omp parallel for reduction(+:histo[:size]) reduction(+:testcount2)
             for (int par = 0; par < NUM_THREADS; par++)
                 for (int i = par + 1; i < n; i += (int) (-log(1. - drand48())*100. * NUM_THREADS), testcount2++)
                     histo[max(0, min(size - 1, (int) floor(array[i * d + j] / epsilon) - (int) first_int))]++;
             long long ref = 0;
             for (int i = 0; i < size; i++)
                 ref += histo[i]*(histo[i] - 1) / 2;
             for (int i = 1; i < size; i++)
                 ref += histo[i - 1] * histo[i];
             printf("%3d %8d %8d %8d %20ld   [%ld", j, (int) first_int, (int) last_int, size, ref, histo[0]);
             costref[j] = ref;
             for (int i = 1; i < size && i < 10; i++)
                 printf(", %ld", histo[i]);
             if (size > 10)
                 printf(", ...] %d %d\n", testcount, testcount2);
             else
                 printf("]\n");
             free(histo);
         }
     }
     int * reorder_rev = (int*) malloc ((d+8)*sizeof(int));
     for (int j=0 ; j<d+8 ; j++)
         reorder_rev[j] = j;
     qsort(reorder_rev, d, sizeof(int), cmp_reorder_dim);
     for (int j=0 ; j<d+8 ; j++)
         printf("%2d %2d %ld\n", j, reorder_rev[j], j<d?costref[reorder_rev[j]]:0);
     reorder_dim = (int*) malloc ((d+8)*sizeof(int));
     for (int j=0 ; j<d+8 ; j++)
         reorder_dim[reorder_rev[j]] = j;
     free(reorder_rev);
     free(costref);
 }

 int main(int argc, char** argv) {

 //    fillperplex(150,4,iris,irissig,5.);
 //    exit(0);
 //    poweriter(150,4,iris,eigenv);
 //    eigenVectorOrdering(150,4,eigenv,iris);
 //    double * proj = createProjection(150,4,eigenv,iris);
 //    fillperplex(150,4,iris,irissig,5.);
 //    for(int i=0 ; i<150 ; i++)
 //        irissig[i] = sqrt(irissig[i]);
 //    int * list = perpjFillList(150,4,proj,irissig,1.);
 //    for(int i=0 ; i<150 ; i++)
 //        printf("%d %f %f %d %d %d\n", i, proj[i], irissig[i], list[i+256], list[i/2+128], list[i/4+64]);
 //    exit(0);
 //      //floydHilloopOutput(40);
 //      testWarshall();
 //    generateSinData(300,3,3,0.05,(double*)malloc(300*3*sizeof(double)));
 //    int n = 2000000; int d = 11; double epsilon = 0.034; char filename[] = "sinclust4-11-5";
 //    int n= 11000000; int d=29; double epsilon = 1.8; char filename[] = "/home/share/higgs/HIGGS_11000000x29.bin";
 //    int n= 11620300; int d=57; double epsilon = 20.; char filename[] = "/home/share/bigcross/bigCross_11620300x57.bin";
 //    int n= 11620300; int d=57; double epsilon = .06 ; char filename[] = "/home/share/bigcross/bigCross_11620300x57_normalized_reorder_0.06.bin";
     int n= 11620300; int d=57; double epsilon = .06 ; char filename[] = "/home/share/bigcross/bigCross_11620300x57_normalized.bin";
 //    int n= 4898400; int d=42; double epsilon = .00001 ; char filename[] = "/home/share/kddcup1999/kddcup1999_4898431x42.bin";

     omp_set_num_threads(NUM_THREADS);
     printf("n=%d d=%d eps=%lf threads=%d omp_threads=%d\n",n,d,epsilon,NUM_THREADS,omp_get_max_threads());
     double * array = (double*) mallocA64((n+7)/8*8 * sizeof (double) * d + 16384);
 //    generateSinData(n,d,5,0.001,array); FILE *ff = fopen("sinclust4-11-5", "w"); fwrite(array, n*d, sizeof(double), ff); fclose(ff);
 //    generateSinDataGeneral(1000,3,1,3,0.001,array); //FILE *ff = fopen("test1000", "w"); fwrite(array, n*d, sizeof(double), ff); fclose(ff);
 //    exit(0);
     FILE * f3 = fopen(filename, "r");
 //    FILE * f3 = fopen("sinclust283", "r");
 //    FILE * f3 = fopen("bmatrix", "r");
     fread(array, n*d, sizeof (double), f3);
     fclose(f3);
 //    vec tvec = _mm512_load_pd(array);
 //    tvec = _mm512_erf_pd(tvec);
 //int testa[]={0,1054,2195,3082,4063,5180,6510,8007,9654,13104,15963,19489,23814,27929,31657,34413,36317,38004,
 //38990,40105,41207,42252,43454,44939,47602,51265,53771,56952,60414,64732,71107,78351,85364,88583,90142,94583,
 //96334,98157,101605,107460,114912,122387,130360,131984,138971,140919,142687,144881,149451,153593,159226,163774,
 //166062,175223,188181,194366,197626,204250,208703,213158,222695,234209,239046,242537,250000};
 //for(int i=0 ; i<64 ; i++) testa[i]=(int)(drand48()*250000);
 //#pragma omp parallel for num_threads(64) proc_bind(close)
 //    for(int i=0 ; i<64 ; i++){
 //        double h=0 ;
 //        for(int j=testa[i]*8*d ; j<testa[i+1]*8*d ; j+=2)
 //            h+=array[j];
 //        printf("%d %d %d %f\n", omp_get_thread_num(), i, testa[i], h);
 //    }

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

 //    printf("vorher Test Hilbert Loop\n");
 //    stop();
 //    stopc();
 //       int i=0 ; int j=0;
 ////       FUR_HILBERT_START(i,j,0,300000,0,300000){
 //       FGF_HILBERT_FOR(i, j, 1 << 30 , 1 << 30, i >=1000000000 && i < 1000000010 && j >= 1000000000 && j < 1000000010, FURHIL_lb0 <=1000000010 && FURHIL_ub0 >=1000000000 && FURHIL_lb1<=1000000010 && FURHIL_ub1 >= 1000000000) {
 ////           if(i >=200000 && i < 200010 && j >= 200000 && j < 200010 || HILLOOP_hilbert % 10000000ULL == 0ULL)
 ////           if(HILLOOP_hilbert == 200000000ULL)
 ////           if(HILLOOP_hilbert % 10000000ULL == 0ULL)
 //           printf("%ld %d %d \n", HILLOOP_hilbert, i,j);
 ////           printf("%ld %d %d %d %d %d %d %d %d\n", HILLOOP_hilbert, i,j, HILLOOP_isize, HILLOOP_jsize, (HILLOOP_J + 1) * (HILLOOP_jcur-HILLOOP_j57), HILLOOP_j57, HILLOOP_J, HILLOOP_jcur);
 ////           if(HILLOOP_hilbert == 0ULL){
 ////               printf("%ld\n", HILLOOP_stop);
 //
 //
 //       }FGF_HILBERT_END(i,j)
 ////       }FUR_HILBERT_END(i,j)
 //    printf("FURHIL | %8.3f | %8.3f\n", stopc(), stop());

     stop();
     stopc();
     outputStatistics(n,d,epsilon,array);
     printf("statistics    | %8.3f | %8.3f\n", stopc(), stop());
 //    transpose_8xd_reorder(n,d,array);
 //    transpose_dx8(n,d,array);
     reorder_dimensions(n,d,array);
     printf("reorder       | %8.3f | %8.3f\n", stopc(), stop());
 //    outputStatistics(n,d,epsilon,array);
 //    printf("statistics    | %8.3f | %8.3f\n", stopc(), stop());
     stop();
     stopc();
     test_ego_loop3(n,d,epsilon,array);
     printf("test_ego_loop3| %8.3f | %8.3f\n", stopc(), stop());
     exit(0);
     test_ego_loop_nested(n,d,epsilon,array);
     printf("test_ego_loop_nested| %8.3f | %8.3f\n", stopc(), stop());
 //x    exit(0);
     long long *counter = corePoints3(n,d,epsilon,array);
     printf("corePoints3| %8.3f | %8.3f\n", stopc(), stop());
     dbscanOnepass(n,d,epsilon,10,array);
     printf("dbscanOnepass| %8.3f | %8.3f\n", stopc(), stop());
     dbscanOnepassUf(n,d,epsilon,10,array);
     printf("dbscanOnepassUf| %8.3f | %8.3f\n", stopc(), stop());
     dbscan(n,d,epsilon,10,array);
     printf("dbscan| %8.3f | %8.3f\n", stopc(), stop());
     dbscan_semaphor(n,d,epsilon,11,array);
     printf("dbscan_semaphor| %8.3f | %8.3f\n", stopc(), stop());
 //    long long * counter = corePoints(n,d,epsilon,array);
 //    printf("corePoints| %8.3f | %8.3f\n", stopc(), stop());
     counter = corePointsMutex(n,d,epsilon,array);
     printf("corePointsMutex| %8.3f | %8.3f\n", stopc(), stop());

     dbscan_simple(n,d,epsilon,10,array);
     printf("dbscanSimple| %8.3f | %8.3f\n", stopc(), stop());


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

 //    dbscan(n,d,epsilon,10,array);
 //    printf("dbscan| %8.3f | %8.3f\n", stopc(), stop());
     exit(0);
 //    long long * counter = corePoints(n,d,epsilon,array);
     for (int i=0 ; i<100 ; i++) printf("%ld ", counter[i]);
     int all=0 ; for (int i=0 ; i<n ; i++)all+=counter[i];
     printf("\nall: %d\n",all);
     printf("Core Points| %8.3f | %8.3f\n", stopc(), stop());
     test_ego_loop(n,d,epsilon,array);
     printf("test loop 1| %8.3f | %8.3f\n", stopc(), stop());
     test_ego_loop(n,d,epsilon,array);
     printf("test loop 1| %8.3f | %8.3f\n", stopc(), stop());
 //    test_ego_loop2(n,d,epsilon,array);
 //    printf("test loop 2| %8.3f | %8.3f\n", stopc(), stop());
 //    long long * counter = corePoints(n,d,epsilon,array);
 //    for (int i=0 ; i<100 ; i++) printf("%ld ", counter[i]);
 //    int all=0 ; for (int i=0 ; i<n ; i++)all+=counter[i];
 //    printf("\nall: %d\n",all);
 //    printf("Core Points| %8.3f | %8.3f\n", stopc(), stop());

     transClosureDens(n,d,epsilon,3,array,counter);
     printf("TransClosur| %8.3f | %8.3f\n", stopc(), stop());
     for(int i=0 ; i<200 ; i++)
         printf("%5d %5d %5d\n", i, listArray[i]->minid, listArray[i]->n) ;
     epsilonGridJoinCanonicAVX512(n, d, epsilon, array);
     printf("EGOJoin 512| %8.3f | %8.3f\n", stopc(), stop());
     double *self = mallocA64(sizeof (double) * n);
     int numStripes = 20;
     int **lower = (int **) malloc (numStripes*sizeof(int*));
     int **upper = (int **) malloc (numStripes*sizeof(int*));
     prepareTwoStripes(n, d, epsilon, array, lower, upper, self);
 //    for (int i=0 ; i<nn/1024 ; i++)
 //        printf("%8d %8d %8d %8d %8d\n", i, lower[0][i+nn/1024], upper[0][i+nn/1024], lower[1][i+nn/1024], upper[1][i+nn/1024]);
     printf("Prepare    | %8.3f | %8.3f\n", stopc(), stop());
     prepareStripes(n, d, numStripes, epsilon, array, lower, upper, self);
 //    for (int i=0 ; i<200 ; i++)
 //        printf("%8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d\n", i, lower[0][i+nn/1024], upper[0][i+nn/1024], lower[1][i+nn/1024],
 //                upper[1][i+nn/1024], lower[2][i+nn/1024], upper[2][i+nn/1024], lower[3][i+nn/1024], upper[3][i+nn/1024], lower[4][i+nn/1024], upper[4][i+nn/1024]);
 //    for (int i=0 ; i<200 ; i++)
 //        printf("%8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d\n", i, lower[0][i+nn/8], upper[0][i+nn/8], lower[1][i+nn/8], 
 //                upper[1][i+nn/8], lower[2][i+nn/8], upper[2][i+nn/8], lower[3][i+nn/8], upper[3][i+nn/8], lower[4][i+nn/8], upper[4][i+nn/8]);
     printf("PreparePar | %8.3f | %8.3f\n", stopc(), stop());
     epsilonGridJoinSimplified(n, d, epsilon, array, numStripes, lower, upper, self);
 //    for(int i=0 ; i<20 ; i++)
 //        printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
 //                testarray[i][0], testarray[i][1], testarray[i][2], testarray[i][3], testarray[i][4], testarray[i][5], testarray[i][6], testarray[i][7], testarray[i][8], testarray[i][9],
 //                testarray[i][10], testarray[i][11], testarray[i][12], testarray[i][13], testarray[i][14], testarray[i][15], testarray[i][16], testarray[i][17], testarray[i][18], testarray[i][19]);
 //    for (int i=0 ; i<25000 ; i++){
 //        for (int s=0 ; s<numStripes ; s++)
 //            for (int j=lower[s][i+nn/8]; j<upper[s][i+nn/8] ; j++)
 //                if(testarray[i][j] != 1)
 //                    printf("ERROR ijs %d %d %d %d lower %d upper %d \n", i, j, s, testarray[i][j], lower[s][i+nn/8], upper[s][i+nn/8]) ;
 //        for (int s=1 ; s<numStripes ; s++)
 //            for (int j=upper[s-1][i+nn/8]; j<lower[s][i+nn/8] ; j++)
 //                if(testarray[i][j] != 0)
 //                    printf("ERROR ijs %d %d %d %d\n", i, j, s, testarray[i][j]) ;
 //    }
     printf("New Simple | %8.3f | %8.3f\n", stopc(), stop());
     JoinCanonicAVX512(n, d, epsilon, array);
     printf("Join AVX512| %8.3f | %8.3f\n", stopc(), stop());
     joinCanonicAVX1(n, d, epsilon, array);
     printf("Join AVX1  | %8.3f | %8.3f\n", stopc(), stop());
     joinCanonic(n, d, epsilon, array);
     return 0;
 }
