cmake_minimum_required(VERSION 3.10)

# set the project name
project(SimISPCMiscTest LANGUAGES C VERSION 1.0)

add_library(SimISPCMultiThreadMisc.o SHARED SimISPCMultiThreadMisc.c)

# dynamic linking: target_link_libraries(SimISPCMiscTest -lpthread)
# static linking:  target_link_libraries(SimISPCMiscTest libpthread.a)

if (TIME_D)
    add_compile_definitions(-DTIME_D=1)
elseif(TIME_THR)
    add_compile_definitions(-DTIME_THR=1)
elseif (TIME_TOT)
   add_compile_definitions(-DTIME_TOT=1)
endif()

# add the executable
add_executable(SimISPCMiscTest SimISPCMiscTest.c)

if (ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES clang-tidy clang-tidy-mp-11)
    if(CLANG_TIDY_EXE)
        message("found clang-tidy")
        set(CLANG_TIDY_COMMAND "${CLANG_TIDY_EXE}" "-checks=-*,modernize-*")
        set(CMAKE_CXX_CLANG_TIDY
          clang-tidy;
          -header-filter="\\.{0,1}.*\\.h$";
          -checks=*;)
        set_target_properties(SimISPCMiscTest PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
    endif()
endif()

target_link_libraries(SimISPCMiscTest SimISPCMultiThreadMisc.o ${stdio_LIBRARY} ${stlib_LIBRARY} ${math__LIBRARY} ${pthread__LIBRARY} ${assert__LIBRARY})

#include_directories(${stdio_INCLUDE_DIRS} ${stlib_INCLUDE_DIRS} ${math_INCLUDE_DIRS} ${pthread_INCLUDE_DIRS} ${assert_INCLUDE_DIRS})
target_include_directories(SimISPCMiscTest INTERFACE ${CMAKE_SOURCE_DIR}/SimISPCMultiThreadMisc.h)

#install(TARGETS SimISPCMiscTest DESTINATION ./)