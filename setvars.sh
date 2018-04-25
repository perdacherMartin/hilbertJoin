export CXX=g++
export CC=gcc

export OMP_NUM_THREADS=64
unset OMP_NESTED
unset OMP_NUM_THREADS
unset OMP_PROC_BIND
unset OMP_PLACES

# export OMP_PLACES=`numactl -H | grep cpus | awk '(NF>3) {for (i = 4; i <= NF; i++) printf "%d,", $i}' | sed 's/.$//'`
# export OMP_NESTED=TRUE
# export OMP_NUM_THREADS=4,64
# export OMP_PROC_BIND=spread,close
#
# # hot teams keep thread pool alive, and remove overhead due to creating/destroying threads
# export KMP_HOT_TEAMS_MAX_LEVEL=2
# export KMP_HOT_TEAMS_MODE=1

# stacksize is usually set too smal on KNL. 16 MB is recommended size on any processor
# export KMP_STACKSIZE=

# bind threads to physical processing units, works only on WINDOWS and Linux not on (OS X*)
# IF KMP_AFFINITY is set, OMP_PROC_BIND and OMP_PLACES will be ignored
# export KMP_AFFINITY=

# controlling time prior a thread goes to sleep, default 200 ms. Setting this to INF can be benefitial for apps which has many parallel regions
# export KMP_BLOCKTIME
