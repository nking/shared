cmake_minimum_required(VERSION 3.10)

# set the project name
project(SimISPCMiscTest LANGUAGES C VERSION 1.0)

add_library(SimISPCMultiThreadMisc.o SHARED SimISPCMultiThreadMisc.c)

# dynamic linking: target_link_libraries(SimISPCMiscTest -lpthread)
# static linking:  target_link_libraries(SimISPCMiscTest libpthread.a)

# add the executable
add_executable(SimISPCMiscTest SimISPCMiscTest.c)

target_link_libraries(SimISPCMiscTest SimISPCMultiThreadMisc.o ${stdio_LIBRARY} ${stlib_LIBRARY} ${math__LIBRARY} ${pthread__LIBRARY} ${assert__LIBRARY})

#include_directories(${stdio_INCLUDE_DIRS} ${stlib_INCLUDE_DIRS} ${math_INCLUDE_DIRS} ${pthread_INCLUDE_DIRS} ${assert_INCLUDE_DIRS})
target_include_directories(SimISPCMiscTest INTERFACE ${CMAKE_SOURCE_DIR}/SimISPCMultiThreadMisc.h>)

install(TARGETS SimISPCMiscTest DESTINATION ./)