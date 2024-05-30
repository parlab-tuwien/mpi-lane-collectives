
cmake_minimum_required(VERSION 3.22)

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")
   
project(multilanecomm
        VERSION 0.1.0
        LANGUAGES C)

find_package(MPI REQUIRED)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")

#set(INCLUDE_PLATFORM_CONFIG_FILE "${CMAKE_SOURCE_DIR}/platform_files/default.cmake"
#    CACHE STRING "Configure project to use a specific platform file")
#include(${INCLUDE_PLATFORM_CONFIG_FILE})

#option(DISABLE_ZERO_COPY "Disable zero-copy collectives [default: enabled]" off)

#SET(MY_COMPILE_FLAGS "")

#if(DISABLE_ZERO_COPY)
#    SET(MY_COMPILE_FLAGS "${MY_COMPILE_FLAGS} -DENABLE_ZEROCOPY=0")
#else()
#    SET(MY_COMPILE_FLAGS "${MY_COMPILE_FLAGS} -DENABLE_ZEROCOPY=1")
#endif()

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_COMPILE_FLAGS}")

set(LANE_COMMON_SRC_FILES
${SRC_DIR}/tuw_lanecoll.c
${SRC_DIR}/tuw_lanecoll.h
${SRC_DIR}/tuw_lane_bcast.c
)

add_executable(lanetest_all
${LANE_COMMON_SRC_FILES}
${SRC_DIR}/tuw_lanecoll_bench_all.c
${SRC_DIR}/tuw_lanecoll_argparse.c
${SRC_DIR}/tuw_lane_utils.c
)
TARGET_LINK_LIBRARIES(lanetest_all MPI::MPI_C)


set(BENCHMARK_NAME ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION})
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${BENCHMARK_NAME}-src")
set(DIST_ARCHIVE "${BENCHMARK_NAME}-dist")

set(CPACK_SOURCE_GENERATOR "TBZ2")
SET(CPACK_SOURCE_IGNORE_FILES
"/CMakeFiles/"
"/_CPack_Packages/"
"/bin/"
"/include/"
"/lib/"
"/scripts/"
"/create_package.sh"
"/Makefile"
"\\\\.tar\\\\.bz2$"
"/\\\\.git/"
"/\\\\.settings/"
"/\\\\..*"
"CMakeCache.txt"
"cmake_install\\\\.cmake"
"CPackConfig\\\\.cmake"
"CPackSourceConfig\\\\.cmake"
"sources\\\\.org"
"${CPACK_SOURCE_IGNORE_FILES}"
)

INCLUDE(CPack)

add_custom_target(pack COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
