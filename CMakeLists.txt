cmake_minimum_required(VERSION 3.6)
project(ego2)

## set gcc compiler, tested with gcc 6.1.0
# export CC=/opt/gcc6.1.0/usr/local/bin/gcc
# export CXX=/opt/gcc6.1.0/usr/local/bin/g++

set(UTIL_SOURCES util/dataIo.cpp util/arguments.cpp util/allocation.cpp  measure/timer.cpp util/chrisutil.cpp)

set(SOURCE_FILES_SELF main.cpp hilbertjoin/egojoin.cpp hilbertjoin/hilloop.cpp ${UTIL_SOURCES})
set(SOURCE_FILES_JOIN mainJoin.cpp hilbertjoin/egojoin.cpp hilbertjoin/hilloop.cpp ${UTIL_SOURCES})

#####################
# build type: Release
#####################

# set(CMAKE_CXX_FLAGS  "-std=c++11 -march=knl -mtune=knl -fpic -ffast-math -O3 -DNDDEBUG -fopenmp")

# adding MKL include directory
include_directories($ENV{MKLROOT}/include)

# set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -march=knl -mtune=knl -ffast-math -fassociative-math -O3 -fopenmp -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lboost_system")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -march=knl -mtune=knl -ffast-math -fassociative-math -O3 -fopenmp -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5")

# xeon-phi
link_directories($ENV{MKLROOT}/lib/intel64)

include_directories($ENV{BOOST_HOME}/include/boost)
link_directories($ENV{BOOST_HOME}/lib)

# mkl libraries
find_library ( mkl_lp64_LIB NAMES libmkl_intel_lp64.a
        PATHS $ENV{MKLROOT} PATH_SUFFIXES lib)
find_library ( mkl_core_LIB NAMES libmkl_core.a
        PATHS $ENV{MKLROOT} PATH_SUFFIXES lib)
find_library ( mkl_thread_LIB NAMES libmkl_intel_thread.a
        PATHS $ENV{MKLROOT} PATH_SUFFIXES lib)
find_library ( mkl_omp_LIB NAMES libiomp5.a
        PATHS $ENV{MKLROOT} PATH_SUFFIXES lib)

####################
# build type: Debug
####################
# cmake -DCMAKE_BUILD_TYPE=Debug ..
# and ignore the warning: "Manually-specified variables were not used by the project: CMAKE_BUID_TYPE"

# set(CMAKE_CXX_FLAGS_DEBUG  "-std=c++11 -march=knl -mtune=knl -fpic -ffast-math -O0 -DNDDEBUG -fopenmp -lboost_system -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5")
set(CMAKE_CXX_FLAGS_DEBUG  "-std=c++11 -march=knl -mtune=knl -fpic -ffast-math -O0 -DNDDEBUG -fopenmp -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5")


# add_executable(egoHilb ${SOURCE_FILES})
add_executable(hilbertJoinCountOnly ${SOURCE_FILES_JOIN})
add_executable(hilbertSelfJoinCountOnly ${SOURCE_FILES_SELF})

target_compile_definitions(hilbertJoinCountOnly PRIVATE -DCOUNT_ONLY)
target_compile_definitions(hilbertSelfJoinCountOnly PRIVATE -DCOUNT_ONLY)

## for a non-verbose version comment out the next two lines
target_compile_definitions(hilbertJoinCountOnly PRIVATE -DOUTPUT)
target_compile_definitions(hilbertSelfJoinCountOnly PRIVATE -DOUTPUT)

if ($ENV{KBLOCK})
    target_compile_definitions(hilbertJoinCountOnly PRIVATE -DKBLOCK=$ENV{KBLOCK})
    target_compile_definitions(hilbertSelfJoinCountOnly PRIVATE -DKBLOCK=$ENV{KBLOCK})
endif()

if ($ENV{OMP_NUM_THREADS})
    target_compile_definitions(hilbertJoinCountOnly PRIVATE -DNUM_THREADS=$ENV{OMP_NUM_THREADS})
    target_compile_definitions(hilbertSelfJoinCountOnly PRIVATE -DNUM_THREADS=$ENV{OMP_NUM_THREADS})
endif()


# message("MKLROOT is $ENV{MKLROOT}")
# message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
# message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
# message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
