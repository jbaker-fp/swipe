# FindFFTW.cmake
find_path(FFTW_INCLUDES fftw3.h)
find_library(FFTW_LIBRARIES NAMES fftw3)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

if(FFTW_FOUND AND NOT TARGET FFTW::FFTW)
  add_library(FFTW::FFTW UNKNOWN IMPORTED)
  set_target_properties(FFTW::FFTW PROPERTIES
    IMPORTED_LOCATION "${FFTW_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDES}")
endif()

mark_as_advanced(FFTW_LIBRARIES FFTW_INCLUDES)
