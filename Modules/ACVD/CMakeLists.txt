#LIB ACVD
MESSAGE(STATUS "Start to parse CMake file for ACVD..")
MESSAGE(STATUS "Current dir: " ${CMAKE_CURRENT_SOURCE_DIR})
MESSAGE(STATUS "Output dir: " ${LIBRARY_OUTPUT_PATH})
					 
INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}/Common
${CMAKE_CURRENT_SOURCE_DIR}/DiscreteRemeshing)
# --------------------------------------------------------------------------
# Library compilation
ADD_LIBRARY(ACVD
${CMAKE_CURRENT_SOURCE_DIR}/Common/RenderWindow.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkDelaunay.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkDualMeshDisplay.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkOFFReader.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkOFFWriter.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkQuadricTools.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkRandomTriangulation.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkSMFReader.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkSurface.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkSurfaceBase.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkVolumeProperties.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkNeighbourhoodComputation.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkCurvatureMeasure.cxx
${CMAKE_CURRENT_SOURCE_DIR}/Common/vtkMyMinimalStandardRandomSequence.cxx
${CMAKE_CURRENT_SOURCE_DIR}/DiscreteRemeshing/vtkDiscreteRemeshing 
${CMAKE_CURRENT_SOURCE_DIR}/DiscreteRemeshing/xmlParser.cpp
${CMAKE_CURRENT_SOURCE_DIR}/DiscreteRemeshing/vtkManifoldSimplification.cxx)
TARGET_LINK_LIBRARIES(ACVD ${VTK_LIBRARIES})

MESSAGE(STATUS "Parse CMake file for ACVD done.")
