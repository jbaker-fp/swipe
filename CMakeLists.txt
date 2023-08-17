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

#--------------------------------------
# Modules
#--------------------------------------
# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# find Eigen
find_package(Eigen3 REQUIRED)
MESSAGE(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")

# find FFTW
find_package(FFTW REQUIRED)
MESSAGE(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")


#--------------------------------------
# Includes
#--------------------------------------
include_directories(
    src
    include
    ${EIGEN3_INCLUDE_DIR}
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
    FFTW::FFTW
    -lblas
    -llapack
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

