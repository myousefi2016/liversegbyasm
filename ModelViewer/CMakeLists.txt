cmake_minimum_required(VERSION 2.8)
PROJECT(ModelViewer)
SET(CMAKE_VERBOSE_MAKEFILE ON)

SET( ModelViewr_SRCS
    main.cpp
	ModelViewer.cpp)

SET( ModelViewr_UIS 
	ModelViewer.ui)

SET( ModelViewr_MOC_HDRS
	ModelViewer.h)
	
QT4_WRAP_UI( ModelViewr_UIS_H ${ModelViewr_UIS} )
QT4_WRAP_CPP( ModelViewr_MOC_SRCS ${ModelViewr_MOC_HDRS} )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# use what QVTK built with
SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")

ADD_EXECUTABLE( ModelViewer
				main.cpp
				${ModelViewr_SRCS} 
				${ModelViewr_UIS_H} 
				${ModelViewr_MOC_SRCS} 
				)
TARGET_LINK_LIBRARIES( ModelViewer
  QVTK
  ${LiverSegSSM_LIRARIES}
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
)

