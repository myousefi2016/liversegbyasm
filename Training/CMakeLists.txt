MESSAGE(STATUS "Start to parse CMake file for training..")
MESSAGE(STATUS "Current dir: " ${CMAKE_CURRENT_SOURCE_DIR})
MESSAGE(STATUS "Output dir: " ${LIBRARY_OUTPUT_PATH})

add_executable (trainingEnginePose trainingEnginePose.cpp)
target_link_libraries (trainingEnginePose  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (trainingEngineShape trainingEngineShape.cpp)
target_link_libraries (trainingEngineShape  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (trainingEngineAppearance trainingEngineAppearance.cpp ${ANN_SRC})
target_link_libraries (trainingEngineAppearance ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (trainingEngineShapeAlignment trainingEngineShapeAlignment.cpp)
target_link_libraries (trainingEngineShapeAlignment  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})


