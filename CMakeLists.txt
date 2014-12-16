CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(LiverSegASM)
SET(CMAKE_VERBOSE_MAKEFILE ON)

MESSAGE(STATUS "Start to parse CMake file for LiverSegSSM..")
MESSAGE(STATUS "Project source dir: " ${PROJECT_SOURCE_DIR})
MESSAGE(STATUS "Project source dir: " ${PROJECT_BINARY_DIR})

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

FIND_PACKAGE(ITK REQUIRED)
include (${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
include (${VTK_USE_FILE})

FIND_PACKAGE(Qt4 REQUIRED)
include (${QT_USE_FILE})

set(LiverSegSSM_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(LiverSegSSM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/Modules/itkAuxiliary
							${PROJECT_SOURCE_DIR}/Modules/kmCommon
							${PROJECT_SOURCE_DIR}/Modules/statismo
							${PROJECT_SOURCE_DIR}/Modules/AdaSegment)

MESSAGE(STATUS "Start to build 3rdParty libraries..")
add_subdirectory(Modules/ACVD)
add_subdirectory(Modules/kmCommon)
add_subdirectory(Modules/AdaSegment)

set(KM_LIBRARIES kmCommon)
include_directories(${LiverSegSSM_INCLUDE_DIR})

add_library(LiverSegmentationAPI STATIC LiverSegmentationAPI.cpp)
target_link_libraries (LiverSegmentationAPI ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES} )

add_executable(LiverSegmentation LiverSegmentation.cpp)
target_link_libraries (LiverSegmentation LiverSegmentationAPI)

#install(TARGETS LiverSegmentationAPI DESTINATION lib EXPORT LiverSegASM-targets)
add_subdirectory(Tests)
add_subdirectory(Training)
add_subdirectory(BuildingModel)
add_subdirectory(Modules/ModelViewer)

MESSAGE(STATUS "Parse CMake file for LiverSegSSM done.")
