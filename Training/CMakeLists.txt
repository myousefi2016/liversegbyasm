CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(TrainingEngine)
SET(CMAKE_VERBOSE_MAKEFILE ON)

add_executable (trainingEnginePose trainingEnginePose.cpp)
target_link_libraries (trainingEnginePose  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (trainingEngineShape trainingEngineShape.cpp)
target_link_libraries (trainingEngineShape  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (trainingEngineAppearance trainingEngineAppearance.cpp ${ANN_SRC})
target_link_libraries (trainingEngineAppearance ${LiverSegSSM_LIRARIES}  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (trainingEngineShapeAlignment trainingEngineShapeAlignment.cpp)
target_link_libraries (trainingEngineShapeAlignment  ${ITK_LIBRARIES} ${VTK_LIBRARIES})


