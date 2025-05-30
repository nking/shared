# mkdir build
# cd build
# cmake ../
# cmake cmake --build .
cmake_minimum_required(VERSION 3.15)

# for tasksys.cpp, need:
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_NAME ISPCMiscTest)
project(${PROJECT_NAME})

add_library(tasksys.o tasksys.cpp)

#find_package(ISPC REQUIRED)
#set(ISPC_SOURCES ispc_function.ispc)
#ispc_compile(ispc_function ${ISPC_SOURCES})

if (TIME_D)
    set(macro_time '-DTIME_D=1')
    add_compile_definitions(-DTIME_D=1)
elseif(TIME_THR)
    set(macro_time '-DTIME_THR=1')
    add_compile_definitions(-DTIME_THR=1)
elseif (TIME_TOT)
    set(macro_time '-DTIME_TOT=1')
   add_compile_definitions(-DTIME_TOT=1)
else()
    ## don't use empty string here:
    set(macro_time )
endif()

## need 32-bit addressing and 8-wide instruction set.
## can find target choices for the individual platform using
## ispc --target --help
## e.g. some target choices that are 8 wide are:
##    sse2-i32x8, sse4-i32x8, avx1-i32x8, avx2-i32x8, neon-i32x8,
##    sse4-i16x8,
##    avx512skx-x8
##
## for code as given in lecture, can use this: --addressing=32 --target=avx2-i32x16
##   but it doesn't exploit use of programIndex and vector lanes for the incremented
##   programIndex (e.g. for N=16 and SIMD width=8, the programIndex=8 round of
##   execution(s) picks up the result from idx=7 and multiplies it to start the
##   multiplication for programIndex=8.

add_custom_command(
    OUTPUT ispc_function.o
    COMMAND ispc -g ${macro_time} -o ispc_function.o --pic -h ${CMAKE_SOURCE_DIR}/ispc_function.h --addressing=32 --target=avx2-i32x8 ${CMAKE_SOURCE_DIR}/ispc_function.ispc
    DEPENDS ispc_function.ispc
)
add_custom_command(
    OUTPUT ispc_function16.o
    COMMAND ispc -g ${macro_time} -o ispc_function16.o -h ${CMAKE_SOURCE_DIR}/ispc_function.h --addressing=32 --target=avx2-i32x16 ${CMAKE_SOURCE_DIR}/ispc_function.ispc
    DEPENDS ispc_function.ispc
)

#message("macro_time=${macro_time}")

add_executable(ISPCMiscTest ISPCMiscTest.c ispc_function.o tasksys.cpp)
add_executable(ISPCMisc16Test ISPCMiscTest.c ispc_function16.o tasksys.cpp)

#target_link_libraries(ISPCMiscTest ispc_function.o)

target_include_directories(ISPCMiscTest INTERFACE ${CMAKE_SOURCE_DIR}/ispc_function.h ${CMAKE_SOURCE_DIR}/../time_log.h)

