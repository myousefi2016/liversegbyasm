cmake_minimum_required(VERSION 2.6)
# This project is designed to be built outside the Insight source tree.
PROJECT(AdaSegment)

# rpath settings
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_INSTALL_PREFIX}/lib/hdf5;")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Setup build locations.
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(HDF5)
IF(HDF5_FOUND)
  MESSAGE(${HDF5_LIBRARIES})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(HDF5_FOUND)

SET(ImageSamplerCode

 ImageSamplers/itkImageRandomSampler.h
 ImageSamplers/itkImageRandomSampler.txx
 ImageSamplers/itkImageRandomSamplerBase.h
 ImageSamplers/itkImageRandomSamplerSparseMask.h
 ImageSamplers/itkImageRandomSamplerSparseMask.txx
 ImageSamplers/itkImageSample.h
 ImageSamplers/itkImageSamplerBase.h
 ImageSamplers/itkImageSamplerBase.txx
 ImageSamplers/itkImageToVectorContainerFilter.h
 ImageSamplers/itkImageToVectorContainerFilter.txx 
 ImageSamplers/itkVectorContainerSource.h
 ImageSamplers/itkVectorContainerSource.txx
 ImageSamplers/itkVectorDataContainer.h
 ImageSamplers/itkVectorDataContainer.txx )
 
 MESSAGE(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/ImageSamplers )

#ADD_EXECUTABLE(AdaLearn AdaLearn.cxx ${ImageSamplerCode})
#ADD_EXECUTABLE(AdaLearn2D AdaLearn2D.cxx ${ImageSamplerCode})
#ADD_EXECUTABLE(Weightedresample NonRandomResample.cxx ${ImageSamplerCode})
#ADD_EXECUTABLE(AdaSegment AdaSegment.cxx)
#ADD_EXECUTABLE(AdaSegment2 AdaSegment2.cxx)
#ADD_EXECUTABLE(AdaSegment2D AdaSegment2D.cxx)
#ADD_EXECUTABLE(ExtractFeatures ExtractFeatures.cpp)
#ADD_EXECUTABLE(AdaLearn2 AdaLearn2.cpp)
#ADD_EXECUTABLE(test test.cxx)
ADD_LIBRARY(AdaSegment STATIC AdaSegmentAPI.cpp)
 
#TARGET_LINK_LIBRARIES(AdaLearn ${ITK_LIBRARIES})
#TARGET_LINK_LIBRARIES(AdaLearn2D ${ITK_LIBRARIES})
#TARGET_LINK_LIBRARIES(Weightedresample ${ITK_LIBRARIES})
#TARGET_LINK_LIBRARIES(AdaSegment ${ITK_LIBRARIES})
#TARGET_LINK_LIBRARIES(AdaSegment2 ${ITK_LIBRARIES})
#TARGET_LINK_LIBRARIES(AdaSegment2D ${ITK_LIBRARIES})
#TARGET_LINK_LIBRARIES(ExtractFeatures ${ITK_LIBRARIES} ${HDF5_LIBRARIES})
#TARGET_LINK_LIBRARIES(AdaLearn2 ${ITK_LIBRARIES})
#TARGET_LINK_LIBRARIES(test ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(AdaSegment ${ITK_LIBRARIES})

include(${PROJECT_SOURCE_DIR}/cmake/GenerateAdaSegmentConfig.cmake)

configure_file(${PROJECT_SOURCE_DIR}/cmake/AdaSegmentConfig.cmake.in
               ${PROJECT_BINARY_DIR}/AdaSegmentConfig.cmake @ONLY IMMEDIATE)

