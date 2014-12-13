CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Tests)
SET(CMAKE_VERBOSE_MAKEFILE ON)

#FIND_PACKAGE(ITK REQUIRED)
#include (${ITK_USE_FILE})

#FIND_PACKAGE(VTK REQUIRED)
#include (${VTK_USE_FILE})

#FIND_PACKAGE(Qt4 REQUIRED)
#include (${QT_USE_FILE})

#FIND_PACKAGE(Statismo REQUIRED)

#include_directories(${LiverSegSSM_ROOT_DIR})
#include_directories(${STATISMO_INCLUDE_DIR})
#link_directories(${STATISMO_LIBRARY_DIR})

#add_executable (testDeformableSimplexMesh3DFilter testDeformableSimplexMesh3DFilter.cpp)
#target_link_libraries (testDeformableSimplexMesh3DFilter  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

#add_executable (testDeformableSimplexMesh3DBalloonForceFilter testDeformableSimplexMesh3DBalloonForceFilter.cpp)
#target_link_libraries (testDeformableSimplexMesh3DBalloonForceFilter  ${ITK_LIBRARIES} )

#add_executable (testGeodesicActiveContourImageFilter testGeodesicActiveContourImageFilter.cpp)
#target_link_libraries (testGeodesicActiveContourImageFilter  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

#add_executable (testMeshToImageElasticRegistration testMeshToImageElasticRegistration.cpp)
#target_link_libraries (testMeshToImageElasticRegistration  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

#add_executable (testImageToImageElasticRegistration testImageToImageElasticRegistration.cpp)
#target_link_libraries (testImageToImageElasticRegistration  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (testTemp testTemp.cpp)
target_link_libraries (testTemp  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (testProfileClustering testProfileClustering.cpp)
target_link_libraries (testProfileClustering  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (testProfileClassifer testProfileClassifer.cpp ${ANN_SRC})
target_link_libraries (testProfileClassifer ${LiverSegSSM_LIRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (testNoCrossDeform testNoCrossDeform.cpp ${ANN_SRC})
target_link_libraries (testNoCrossDeform ${LiverSegSSM_LIRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (testSSM testSSM.cpp ${ANN_SRC})
target_link_libraries (testSSM ${LiverSegSSM_LIRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (testProfileContainer testProfileContainer.cpp ${ANN_SRC})
target_link_libraries (testProfileContainer ${LiverSegSSM_LIRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})