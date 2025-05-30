# mkdir build
# cd build
# cmake ../
# cmake cmake --build .
cmake_minimum_required(VERSION 3.15)

#set(PROJECT_NAME SIMDMiscTest)
project(SIMDMiscTest LANGUAGES C VERSION 1.0)

include_directories(${stdio_INCLUDE_DIRS} ${math_INCLUDE_DIRS} ${pthread_INCLUDE_DIRS} ${assert_INCLUDE_DIRS} ${immintrin_INCLUDE_DIRS})

#-march=native
add_library(SIMDMisc.o SHARED SIMDMisc.c)
target_compile_options(SIMDMisc.o PRIVATE -mavx2 -mfma)

#message("CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")
#message("CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")

## enable ONE of these 3 at any time to make timing logs
#message("TIME_TOT=${TIME_TOT}")

if (TIME_D)
    add_compile_definitions(-DTIME_D=1)
elseif(TIME_THR)
    add_compile_definitions(-DTIME_THR=1)
elseif (TIME_TOT)
   add_compile_definitions(-DTIME_TOT=1)
endif()

add_executable(SIMDMiscTest SIMDMiscTest.c)

target_compile_options(SIMDMiscTest PRIVATE -mavx2 -mfma)

target_link_libraries(SIMDMiscTest SIMDMisc.o ${stdio_LIBRARY} ${immintrin_LIBRARY} ${math__LIBRARY} ${pthread__LIBRARY} ${assert__LIBRARY})

#message("CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}")
target_include_directories(SIMDMiscTest INTERFACE ${CMAKE_SOURCE_DIR}/SIMDMisc.h ${CMAKE_SOURCE_DIR}/../time_log.h)

if (ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES clang-tidy clang-tidy-mp-11)
    if(CLANG_TIDY_EXE)
        message("found clang-tidy")
        set(CLANG_TIDY_COMMAND "${CLANG_TIDY_EXE}" "-checks=-*,modernize-*")
        set(CMAKE_CXX_CLANG_TIDY
          clang-tidy;
          -header-filter="\\.{0,1}.*\\.h$";
          -checks=*;)
        set_target_properties(SIMDMiscTest PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
    endif()
endif()

