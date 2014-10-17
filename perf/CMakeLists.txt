# CMake entry point
cmake_minimum_required (VERSION 2.6)

project (nbody)

# Specify bin path
set (EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

find_package (OpenGL REQUIRED)

# Compile external dependencies 
add_subdirectory (libs/glfw-3.0.4)

# Compile glew library
add_library( 
	glew STATIC
	libs/glew-1.11.0/src/glew.c
	libs/glew-1.11.0/include
)

# Add include from external dependencies
include_directories (libs/glfw-3.0.4/include)
include_directories (libs/glm-0.9.5.4)
include_directories (libs/glew-1.11.0/include)

# Génération de la liste des fichiers sources
file (GLOB_RECURSE source_files src/*)

# Déclaration de l'exécutable
add_executable (nbody ${source_files})

# Add librairies for linking
target_link_libraries(nbody glew glfw ${GLFW_LIBRARIES})

# Activate C++11
if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
endif()
