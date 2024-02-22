cmake_minimum_required(VERSION 3.16)
project(swipe)
include(ExternalProject)
enable_language(Fortran)

#--------------------------------------
# C++ standard
#--------------------------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#--------------------------------------
# Compiler Flags
#--------------------------------------
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CUSTOM_COMPILER_FLAGS}")
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

#--------------------------------------
# Dependencies
#--------------------------------------

# Sets options for pffft
set(PFFFT_USE_TYPE_FLOAT OFF)
set(PFFFT_USE_TYPE_DOUBLE ON)
set(INSTALL_PFFFT OFF)
set(PFFFT_USE_BENCH_FFTW OFF)
set(PFFFT_USE_BENCH_GREEN OFF)
set(PFFFT_USE_BENCH_KISS OFF)
set(PFFFT_USE_BENCH_POCKET OFF)
set(PFFFT_USE_BENCH_MKL OFF)
set(PFFFT_USE_FFTPACK OFF)

# pffft dependency
add_subdirectory(pffft)

# openblas dependency for local builds and cloud deployments
if (IS_LOCAL_BUILD OR IS_CLOUD_DEPLOYMENT)
    set(OPENBLAS_INSTALL_DIR "${CMAKE_BINARY_DIR}/external/openblas")
    ExternalProject_Add(
        OpenBLAS
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/openblas"
        PREFIX "${CMAKE_BINARY_DIR}/external"
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${OPENBLAS_INSTALL_DIR} -DBUILD_STATIC_LIBS=ON
        BUILD_COMMAND make
        INSTALL_COMMAND make install
    )
    set(OpenBLAS_LIB "${OPENBLAS_INSTALL_DIR}/lib/libopenblas.a")
endif()

#--------------------------------------
# Includes
#--------------------------------------
include_directories(
    src
    include
)

if (NOT EIGEN3_INCLUDE_DIR)
    set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen")
endif()

#--------------------------------------
# Swipe Library
#--------------------------------------
add_library(swipe STATIC
    src/swipe.cpp
    src/fitpack.f
    )
target_include_directories(swipe PRIVATE pfft ${EIGEN3_INCLUDE_DIR})
target_compile_options(swipe PRIVATE -O3 -fPIC)

if (IS_LOCAL_BUILD OR IS_CLOUD_DEPLOYMENT)
    # use openblas in local builds and cloud deployments
    add_dependencies(swipe OpenBLAS)
    target_link_libraries(swipe PRIVATE
        -lgfortran
        ${OpenBLAS_LIB}
        PFFFT
        )
else()
    # use system blas and lapack in firmware builds
    target_link_libraries(swipe PRIVATE
        -lgfortran
        -lblas
        -llapack
        PFFFT
        )
endif()

if(BUILD_BINARY)
    #--------------------------------------
    # Binary
    #--------------------------------------
    add_executable(Swipe_Main src/main.cpp)
    target_link_libraries(Swipe_Main PRIVATE swipe)
    target_compile_options(Swipe_Main PRIVATE -O3)
    set_target_properties(Swipe_Main PROPERTIES OUTPUT_NAME swipe)
endif()

