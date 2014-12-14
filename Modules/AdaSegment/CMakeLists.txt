#LIB AdaSegment
MESSAGE(STATUS "Start to parse CMake file for AdaSegment..")
MESSAGE(STATUS "Current dir: " ${CMAKE_CURRENT_SOURCE_DIR})
MESSAGE(STATUS "Output dir: " ${LIBRARY_OUTPUT_PATH})
	 
INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_CURRENT_SOURCE_DIR}/ImageSamplers)

add_executable(AdaLearn2 
${CMAKE_CURRENT_SOURCE_DIR}/AdaLearn2.cxx)
TARGET_LINK_LIBRARIES(AdaLearn2 ${ITK_LIBRARIES})

add_executable(AdaSegment2 
${CMAKE_CURRENT_SOURCE_DIR}/AdaSegment2.cxx)
TARGET_LINK_LIBRARIES(AdaSegment2 ${ITK_LIBRARIES})

MESSAGE(STATUS "Parse CMake file for AdaSegment done.")

