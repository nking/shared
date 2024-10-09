# mkdir build
# cd build
# cmake ../
# cmake cmake --build .
cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME ISPCMiscTest)
project(${PROJECT_NAME})

#find_package(ISPC REQUIRED)
#set(ISPC_SOURCES ispc_function.ispc)
#ispc_compile(ispc_function ${ISPC_SOURCES})

add_custom_command(
    OUTPUT ispc_function.o
    COMMAND ispc -o ispc_function.o -h ${CMAKE_SOURCE_DIR}/ispc_function.h ${CMAKE_SOURCE_DIR}/ispc_function.ispc
    DEPENDS ispc_function.ispc
)

add_executable(ISPCMiscTest ISPCMiscTest.c ispc_function.o)

#target_link_libraries(ISPCMiscTest ispc_function.o)

target_include_directories(ISPCMiscTest INTERFACE ${CMAKE_SOURCE_DIR}/ispc_function.h>)

