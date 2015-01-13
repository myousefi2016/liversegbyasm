MESSAGE(STATUS "Start to parse CMake file for testing..")
MESSAGE(STATUS "Current dir: " ${CMAKE_CURRENT_SOURCE_DIR})
MESSAGE(STATUS "Output dir: " ${LIBRARY_OUTPUT_PATH})

add_library(testLib testLib.cpp)
target_link_libraries (testLib ${ITK_LIBRARIES} ${KM_LIBRARIES})

add_executable (testTemp testTemp.cpp)
target_link_libraries (testTemp  testLib ${KM_LIBRARIES})

add_executable (testProfileClustering testProfileClustering.cpp)
target_link_libraries (testProfileClustering  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (testProfileClassifer testProfileClassifer.cpp)
target_link_libraries (testProfileClassifer ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (testNoCrossDeform testNoCrossDeform.cpp)
target_link_libraries (testNoCrossDeform ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (testSSM testSSM.cpp)
target_link_libraries (testSSM ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (testProfileContainer testProfileContainer.cpp)
target_link_libraries (testProfileContainer ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (testDeformableSimplexMesh3DWithShapePriorFilter testDeformableSimplexMesh3DWithShapePriorFilter.cpp)
target_link_libraries (testDeformableSimplexMesh3DWithShapePriorFilter ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (testACVD testACVD.cpp)
target_link_libraries (testACVD ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})

add_executable (testGlobalConfig testGlobalConfig.cpp)
target_link_libraries (testGlobalConfig ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${KM_LIBRARIES})