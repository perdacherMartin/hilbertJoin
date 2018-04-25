cmake_minimum_required(VERSION 3.6)
project(ego2)

## set gcc compiler, tested with gcc 6.1.0
# export CC=/opt/gcc6.1.0/usr/local/bin/gcc
# export CXX=/opt/gcc6.1.0/usr/local/bin/g++

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp hilbertjoin/egojoin.cpp util/arguments.cpp util/dataIo.cpp util/chrisutil.cpp hilbertjoin/hilloop.cpp measure/energy.cpp measure/papicalls.cpp measure/timer.cpp util/allocation.cpp)

#####################
# build type: Release
#####################

# set(CMAKE_CXX_FLAGS  "-std=c++11 -march=knl -mtune=knl -fpic -ffast-math -O3 -DNDDEBUG -fopenmp")
set(CMAKE_CXX_FLAGS  "-std=c++11 -march=knl -mtune=knl -ffast-math -fassociative-math -O3 -DNDDEBUG -fopenmp -lmemkind")


####################
# build type: Debug
####################
# cmake -DCMAKE_BUILD_TYPE=Debug ..
# and ignore the warning: "Manually-specified variables were not used by the project: CMAKE_BUID_TYPE"
set(CMAKE_CXX_FLAGS_DEBUG  "-std=c++11 -march=knl -mtune=knl -fpic -ffast-math -O0 -DNDDEBUG -fopenmp")

# papi library
set(PAPI_PREFIX "/usr/local")
find_library(PAPI_LIBRARIES
    # Pick the static library first for easier run-time linking.
    NAMES libpapi.a papi PATHS ${PAPI_PREFIX}/lib
)


add_executable(egoHilb ${SOURCE_FILES})
add_definitions(-DCANOLOOP)
target_link_libraries(egoHilb ${PAPI_LIBRARIES})

add_executable(egoCano ${SOURCE_FILES})
target_link_libraries(egoCano ${PAPI_LIBRARIES})


# message("MKLROOT is $ENV{MKLROOT}")
# message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
# message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
# message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
