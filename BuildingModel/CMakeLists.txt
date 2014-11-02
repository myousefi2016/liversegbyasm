CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(BuildingModel)
SET(CMAKE_VERBOSE_MAKEFILE ON)

add_executable (BuildStatisticalShapeModel BuildStatisticalShapeModel.cpp)
target_link_libraries (BuildStatisticalShapeModel ${LiverSegSSM_LIRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (BuildStatisticalAppearanceModel BuildStatisticalAppearanceModel.cpp)
target_link_libraries (BuildStatisticalAppearanceModel  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (BuildStatisticalPoseModel BuildStatisticalPoseModel.cpp)
target_link_libraries (BuildStatisticalPoseModel  ${ITK_LIBRARIES} ${VTK_LIBRARIES})


