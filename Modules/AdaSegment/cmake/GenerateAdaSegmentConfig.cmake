# Generate the AdaSegmentConfig.cmake file in the build tree.
# The file tells external projects how to use AdaSegment.

# Settings specific to build trees
#
#
set(AdaSegment_USE_FILE_CONFIG ${AdaSegment_BINARY_DIR}/UseAdaSegment.cmake)
set(AdaSegment_INCLUDE_DIRS_CONFIG
  ${AdaSegment_SOURCE_DIR}
  ${AdaSegment_SOURCE_DIR}/ImageSamplers
  )
set(AdaSegment_LIBRARY_DIRS_CONFIG
  ${AdaSegment_BINARY_DIR}
  )
set(AdaSegment_BINARY_DIR_CONFIG
  ${AdaSegment_BINARY_DIR}
  )
set(AdaSegment_LIBRARIES_CONFIG
  AdaSegment
  )

configure_file(
  ${AdaSegment_SOURCE_DIR}/cmake/AdaSegmentConfig.cmake.in
  ${AdaSegment_BINARY_DIR}/AdaSegmentConfig.cmake
  @ONLY IMMEDIATE
  )

configure_file(
  ${AdaSegment_SOURCE_DIR}/cmake/UseAdaSegment.cmake.in
  ${AdaSegment_USE_FILE_CONFIG}
  @ONLY IMMEDIATE
  )

