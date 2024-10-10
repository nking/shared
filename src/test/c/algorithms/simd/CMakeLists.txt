# mkdir build
# cd build
# cmake ../
# cmake cmake --build .
cmake_minimum_required(VERSION 3.15)

#set(PROJECT_NAME SIMDMiscTest)
project(SIMDMiscTest LANGUAGES C VERSION 1.0)

include_directories(${stdio_INCLUDE_DIRS} ${math_INCLUDE_DIRS} ${pthread_INCLUDE_DIRS} ${assert_INCLUDE_DIRS} ${immintrin_INCLUDE_DIRS})

#-march=native
add_executable(SIMDMiscTest SIMDMiscTest.c SIMDMisc.c)

target_compile_options(SIMDMiscTest PRIVATE -mavx2)

target_link_libraries(SIMDMiscTest ${stdio_LIBRARY} ${immintrin_LIBRARY} ${math__LIBRARY} ${pthread__LIBRARY} ${assert__LIBRARY})

target_include_directories(SIMDMiscTest INTERFACE ${CMAKE_SOURCE_DIR}/SIMDMisc.h>)


