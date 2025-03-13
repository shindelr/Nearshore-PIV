#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "png_shared" for configuration "Release"
set_property(TARGET png_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(png_shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libpng16.16.47.0.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libpng16.16.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS png_shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_png_shared "${_IMPORT_PREFIX}/lib/libpng16.16.47.0.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
