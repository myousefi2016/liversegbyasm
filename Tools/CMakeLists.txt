CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Tools)
SET(CMAKE_VERBOSE_MAKEFILE ON)

add_executable (toolFillBinarySliceHole toolFillBinarySliceHole.cpp)
target_link_libraries (toolFillBinarySliceHole  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (toolMeshStatistics toolMeshStatistics.cpp)
target_link_libraries (toolMeshStatistics  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (toolTriangleMeshToSimplexMeshConverter toolTriangleMeshToSimplexMeshConverter.cpp)
target_link_libraries (toolTriangleMeshToSimplexMeshConverter  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (toolSimplexMeshAdaptTopology toolSimplexMeshAdaptTopology.cpp)
target_link_libraries (toolSimplexMeshAdaptTopology  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (toolSimplexMeshToTriangleMeshConverter toolSimplexMeshToTriangleMeshConverter.cpp)
target_link_libraries (toolSimplexMeshToTriangleMeshConverter  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (toolDownsamplePolydata toolDownsamplePolydata.cpp)
target_link_libraries (toolDownsamplePolydata  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (toolComputeMeshGeometry toolComputeMeshGeometry.cpp)
target_link_libraries (toolComputeMeshGeometry  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable (toolProcessForAdboost toolProcessForAdboost.cpp)
target_link_libraries (toolProcessForAdboost  ${ITK_LIBRARIES} ${VTK_LIBRARIES})