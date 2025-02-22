#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "blas" for configuration "RELEASE"
set_property(TARGET blas APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(blas PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libblas.so.3.7.0"
  IMPORTED_SONAME_RELEASE "libblas.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS blas )
list(APPEND _IMPORT_CHECK_FILES_FOR_blas "${_IMPORT_PREFIX}/lib64/libblas.so.3.7.0" )

# Import target "lapack" for configuration "RELEASE"
set_property(TARGET lapack APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(lapack PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "blas"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/liblapack.so.3.7.0"
  IMPORTED_SONAME_RELEASE "liblapack.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS lapack )
list(APPEND _IMPORT_CHECK_FILES_FOR_lapack "${_IMPORT_PREFIX}/lib64/liblapack.so.3.7.0" )

# Import target "tmglib" for configuration "RELEASE"
set_property(TARGET tmglib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(tmglib PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "lapack"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libtmglib.so"
  IMPORTED_SONAME_RELEASE "libtmglib.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS tmglib )
list(APPEND _IMPORT_CHECK_FILES_FOR_tmglib "${_IMPORT_PREFIX}/lib64/libtmglib.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
