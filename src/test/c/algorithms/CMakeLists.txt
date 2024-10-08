cmake_minimum_required(VERSION 3.10)

# set the project name
project(SimISPCMiscTest LANGUAGES C VERSION 1.0)

#add_library(SimISPCMiscTest SimISPCMiscTest.c)

# or dynamic linking: target_link_libraries(SimISPCMiscTest -lpthread)
# or static linking:  target_link_libraries(SimISPCMiscTest libpthread.a)

#find_package(stdio REQUIRED)
#find_package(stdlib REQUIRED)
#find_package(math REQUIRED)
#find_package(pthread REQUIRED)
#find_package(assert REQUIRED)
include_directories(${stdio_INCLUDE_DIRS} ${stlib_INCLUDE_DIRS} ${math_INCLUDE_DIRS} ${pthread_INCLUDE_DIRS} ${assert_INCLUDE_DIRS})

# add the executable
add_executable(SimISPCMiscTest SimISPCMiscTest.c)

target_link_libraries(SimISPCMiscTest ${stdio_LIBRARY} ${stlib_LIBRARY} ${math__LIBRARY} ${pthread__LIBRARY} ${assert__LIBRARY})

# or use ${CMAKE_INSTALL_RUNSTATEDIR}
message("CMAKE_INSTALL_LOCALSTATEDIR=${CMAKE_INSTALL_LOCALSTATEDIR}")
message("CMAKE_INSTALL_RUNSTATEDIR=${CMAKE_INSTALL_RUNSTATEDIR}")
message("CMAKE_INSTALL_RUNSTATEDIR=${CMAKE_INSTALL_RUNSTATEDIR}")

install(TARGETS SimISPCMiscTest DESTINATION ./)