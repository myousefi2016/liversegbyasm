CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(LiverSegASM)
SET(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

option(BUILD_CUDA_LIB "Build CUDA library" OFF)
option(ITK_USE_GPU "GPU acceleration via OpenCL" ON)

FIND_PACKAGE(ITK REQUIRED)
include (${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
include (${VTK_USE_FILE})

FIND_PACKAGE(Qt4 REQUIRED)
include (${QT_USE_FILE})

#FIND_PACKAGE(HDF5 REQUIRED)
#message(${HDF5_INCLUDE_DIR})
#message(${HDF5_LIBRARIES})
#include_directories(${HDF5_INCLUDE_DIR})

FIND_PACKAGE(flann REQUIRED)
include (${FLANN_USE_FILE})

#set(CUDA_LIBRARIES "")
#set(OPENCL_LIBRARIES "")

if (BUILD_CUDA_LIB)
    find_package(CUDA)
    if (CUDA_FOUND)
        message(STATUS "CUDA found (include: ${CUDA_INCLUDE_DIRS}, lib: ${CUDA_LIBRARIES})")
        include_directories(${CUDA_INCLUDE_DIRS})
		set(CUDA_LIBRARIES ${CUDA_LIBRARIES})
    else(CUDA_FOUND)
        message(STATUS "CUDA not found, CUDA library will not be built")
        set(BUILD_CUDA_LIB OFF)
    endif(CUDA_FOUND)
	
	FIND_PACKAGE(OpenCL REQUIRED)
	message(${OPENCL_INCLUDE_DIRS})
	message(${OPENCL_LIBRARIES})
	include_directories(${OPENCL_INCLUDE_DIRS})
endif(BUILD_CUDA_LIB)

set(LiverSegSSM_ROOT_DIR ${LiverSegASM_SOURCE_DIR})
set(LiverSegSSM_INCLUDE_DIR ${LiverSegASM_BINARY_DIR}
							${LiverSegASM_SOURCE_DIR}/include
							${LiverSegASM_SOURCE_DIR}/include/itk
							${LiverSegASM_SOURCE_DIR}/include/km
							${LiverSegASM_SOURCE_DIR}/include/kmstatismo
							${LiverSegASM_SOURCE_DIR}/include/AdaSegment
							${LiverSegASM_SOURCE_DIR}/include/AdaSegment/ImageSamplers
							${OPENCL_INCLUDE_DIRS})
							
include_directories(${LiverSegSSM_INCLUDE_DIR})
link_directories( ${LiverSegSSM_LIBRARY_DIR})

link_libraries(${ITK_LIBRARIES} ${VTK_LIBRARIES} ${CUDA_LIBRARIES} ${OPENCL_LIBRARIES} ${FLANN_LIBRARIES} ${HDF5_LIBRARIES})

#include (GenerateExportHeader)
#add_compiler_export_flags()
add_library(LiverSegmentationAPI STATIC LiverSegmentationAPI.cpp)
#GENERATE_EXPORT_HEADER( LiverSegmentationAPI
#             BASE_NAME LiverSegmentationAPI
#             EXPORT_MACRO_NAME LiverSegmentationAPI_EXPORT
#             EXPORT_FILE_NAME LiverSegmentationAPI_Export.h
#             STATIC_DEFINE LiverSegmentationAPI_BUILT_AS_STATIC
#)

target_link_libraries (LiverSegmentationAPI)

install(TARGETS LiverSegmentationAPI DESTINATION lib EXPORT LiverSegASM-targets)


add_executable(LiverSegmentation LiverSegmentation.cpp)
target_link_libraries (LiverSegmentation LiverSegmentationAPI)

add_executable(TestAdaboostClassifer TestAdaboostClassifer.cpp)
target_link_libraries (TestAdaboostClassifer)

add_executable(ProcessForAdboost ProcessForAdboost.cpp)
target_link_libraries (ProcessForAdboost)

add_executable(AdaLearn2 include/AdaSegment/AdaLearn2.cpp)
target_link_libraries (AdaLearn2)

add_executable(AdaSegment2 include/AdaSegment/AdaSegment2.cxx)
target_link_libraries (AdaSegment2)

#add_executable(TempProcess TempProcess.cpp)
#target_link_libraries (TempProcess)

add_subdirectory(Tools)
#add_subdirectory(Tests)
add_subdirectory(Training)
add_subdirectory(BuildingModel)
add_subdirectory(ModelViewer)

#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lib DESTINATION ${CMAKE_BINARY_DIR})
