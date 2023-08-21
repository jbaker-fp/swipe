cmake_minimum_required(VERSION 3.16)
project(swipe)
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

# # find Eigen
# find_package(Eigen3 REQUIRED)
# MESSAGE(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")

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

# Disable uninstall target for this submodule
add_subdirectory(pffft)


#--------------------------------------
# Includes
#--------------------------------------
include_directories(
    src
    include
    ${EIGEN_INCLUDE_DIRS}
)


#--------------------------------------
# Swipe Library
#--------------------------------------
add_library(swipe STATIC 
    src/swipe.cpp 
    src/fitpack.f
    )
target_include_directories(swipe PUBLIC include)
target_link_libraries(swipe PRIVATE 
    -lblas
    -llapack
    PFFFT
    )
target_compile_options(swipe PRIVATE -O3 -fPIC)

if(BUILD_BINARY)
    #--------------------------------------
    # Binary
    #--------------------------------------
    add_executable(Swipe_Main src/main.cpp)
    target_link_libraries(Swipe_Main PRIVATE swipe)
    target_compile_options(Swipe_Main PRIVATE -O3)
    set_target_properties(Swipe_Main PROPERTIES OUTPUT_NAME swipe)
endif()

