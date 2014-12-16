#LIB kmCommon
MESSAGE(STATUS "Start to parse CMake file for kmCommon..")
MESSAGE(STATUS "Current dir: " ${CMAKE_CURRENT_SOURCE_DIR})
MESSAGE(STATUS "Output dir: " ${LIBRARY_OUTPUT_PATH})

MESSAGE(STATUS ${LiverSegSSM_SOURCE_DIR}/Modules/itkAuxiliary)
include_directories(${LiverSegSSM_SOURCE_DIR}/Modules/itkAuxiliary)

add_library(kmCommon
${CMAKE_CURRENT_SOURCE_DIR}/kmUtility.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmVtkItkUtility.cpp)
TARGET_LINK_LIBRARIES(kmCommon ${ITK_LIBRARIES} ${VTK_LIBRARIES})

MESSAGE(STATUS "Parse CMake file for kmCommon done.")

