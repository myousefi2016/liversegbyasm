MESSAGE(STATUS "Start to parse CMake file for build models..")
MESSAGE(STATUS "Current dir: " ${CMAKE_CURRENT_SOURCE_DIR})
MESSAGE(STATUS "Output dir: " ${LIBRARY_OUTPUT_PATH})

add_executable (BuildStatisticalShapeModel BuildStatisticalShapeModel.cpp)
target_link_libraries (BuildStatisticalShapeModel ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (BuildStatisticalAppearanceModel BuildStatisticalAppearanceModel.cpp)
target_link_libraries (BuildStatisticalAppearanceModel  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (BuildStatisticalPoseModel BuildStatisticalPoseModel.cpp)
target_link_libraries (BuildStatisticalPoseModel  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})


