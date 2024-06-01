cmake_minimum_required(VERSION 3.22)

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")
   
project(multilanecomm
        VERSION 0.1.0
        LANGUAGES C)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

find_package(MPI REQUIRED)

#set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
#set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")

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

# Configure the installation paths
set(LANEC_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LANEC_INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR})
set(LANEC_INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR})
set(LANEC_INSTALL_LIBEXECDIR ${CMAKE_INSTALL_LIBEXECDIR})
set(LANEC_INSTALL_DATADIR ${CMAKE_INSTALL_DATADIR})
set(LANEC_INSTALL_CMAKEDIR ${CMAKE_INSTALL_DATADIR}/cmake/lanec)

STRING(COMPARE EQUAL ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} LANEC_BUILD_INPLACE)
if(NOT LANEC_BUILD_INPLACE)
    include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}")
    include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/include")
    include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/src")
endif(NOT LANEC_BUILD_INPLACE)
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}/include")
include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}/src")

include_directories("${CMAKE_SOURCE_DIR}/include")

add_library(lanec SHARED
        include/tuw_lanecoll.h
        src/tuw_lanecoll.c
        src/tuw_lane_bcast.c
        src/tuw_lanecoll_internal.h
)
add_library(LANEC::lanec ALIAS lanec)
set_target_properties(lanec PROPERTIES PUBLIC_HEADER include/tuw_lanecoll.h)

add_executable(lanetest_all
        examples/tuw_lanecoll_bench_all.c
)
TARGET_LINK_LIBRARIES(lanetest_all lanec MPI::MPI_C)

# make sure the include directory is part of the interface
target_include_directories(lanec
        INTERFACE
        $<INSTALL_INTERFACE:${LANEC_INSTALL_INCLUDEDIR}>
        $<$<NOT:${LANEC_BUILD_INPLACE}>:
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)

target_link_libraries(lanec PUBLIC MPI::MPI_C)

#set_target_properties(mpi_schedule_lib PROPERTIES
#        PUBLIC_HEADER "include/tuw_schedule_collectives.h")

configure_package_config_file(lanec-config.cmake.in
        ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/lanec-config.cmake
        INSTALL_DESTINATION ${LANEC_INSTALL_CMAKEDIR}
        PATH_VARS LANEC_INSTALL_LIBDIR
        LANEC_INSTALL_INCLUDEDIR
        LANEC_INSTALL_BINDIR
        LANEC_INSTALL_LIBEXECDIR
        LANEC_INSTALL_CMAKEDIR)

install(TARGETS lanec
        EXPORT lanec-targets
        LIBRARY DESTINATION ${LANEC_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${LANEC_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${LANEC_INSTALL_BINDIR}
        PUBLIC_HEADER DESTINATION ${LANEC_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_SOURCE_DIR}${CMAKE_FILES_DIRECTORY}/lanec-config.cmake
        DESTINATION ${LANEC_INSTALL_CMAKEDIR})

install(EXPORT lanec-targets # rules to export the targets built here
        FILE LANECTargets.cmake
        NAMESPACE LANEC::
        DESTINATION ${LANEC_INSTALL_CMAKEDIR})

