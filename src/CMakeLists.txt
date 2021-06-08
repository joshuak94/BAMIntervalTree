cmake_minimum_required (VERSION 3.8)

# An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_lib" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE seqan3::seqan3)
# target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE libhts.a)
# target_link_libraries ("${PROJECT_NAME}_lib" INTERFACE lzma)
target_include_directories ("${PROJECT_NAME}_lib" INTERFACE ../include)

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" "${PROJECT_NAME}_lib")
