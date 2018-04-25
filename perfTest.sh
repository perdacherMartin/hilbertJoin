THREADS=(8 16 24 32 40 48 56 64 96 128 160 192 224 256)
# THREADS=(96 128 160 192 224 256)

REPETITIONS=5
BLOCKSIZE=5000

RESULTPATH="/home/martin/projects/EGO/join01/results"

for epsilon in 0.1 0.2 0.3
do
    export RESULTFILE="$RESULTPATH/egoHilb_$epsilon.csv"
    echo "N;D;THREADS;EPS;TIME;COUNTS" >$RESULTFILE

    for t in ${THREADS[@]}
    do
      echo "$t"
      export OMP_NUM_THREADS=$t
      for i in $(seq 1 $REPETITIONS)
      do
         ./build/egoHilb -n 200000 -e $epsilon -d 64 -t $t -f /home/share/test/bmatrix_x_200000x64.bin -b >>$RESULTFILE
    #     numactl --membind=1,2,3,4,5,6,7,8
      done
    done
done

cat $RESULTFILE | mail -s "egoHilbert done!" -r "root@ivanhoe.dm.univie.ac.at" "martin.perdacher@univie.ac.at"

RESULTPATH="/home/martin/projects/EGO/join01/results"

for epsilon in 0.1 0.2 0.3
do
    export RESULTFILE="$RESULTPATH/egoCano_$epsilon.csv"
    echo "N;D;THREADS;EPS;TIME;COUNTS" >$RESULTFILE

    for t in ${THREADS[@]}
    do
      echo "$t"
      export OMP_NUM_THREADS=$t
      for i in $(seq 1 $REPETITIONS)
      do
         ./build/egoCano -n 200000 -e $epsilon -d 64 -t $t -f /home/share/test/bmatrix_x_200000x64.bin -b >>$RESULTFILE
    #     numactl --membind=1,2,3,4,5,6,7,8
      done
    done
done

cat $RESULTFILE | mail -s "egoCano done!" -r "root@ivanhoe.dm.univie.ac.at" "martin.perdacher@univie.ac.at"
