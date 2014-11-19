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

#add_executable (TestDeformSeg TestDeformSeg.cpp)
#target_link_libraries (TestDeformSeg  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

#add_executable (DeformableModel1 DeformableModel1.cpp)
#target_link_libraries (DeformableModel1  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

#add_executable (TestKNNClassifer TestKNNClassifer.cpp ${ANN_SRC})
#target_link_libraries (TestKNNClassifer ${LiverSegSSM_LIRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (TestAdaboostClassifer TestAdaboostClassifer.cpp ${ANN_SRC})
target_link_libraries (TestAdaboostClassifer ${LiverSegSSM_LIRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (TestComposedFitting TestComposedFitting.cpp ${ANN_SRC})
target_link_libraries (TestComposedFitting ${LiverSegSSM_LIRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})