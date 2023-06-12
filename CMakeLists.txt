# CMakeList.txt : CMake project for GLU, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.25 FATAL_ERROR)

string(TIMESTAMP TODAY "%Y%m%d")
string(TIMESTAMP YEAR "%Y")

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project (GLU VERSION ${TODAY} LANGUAGES C CXX)

include(${PROJECT_SOURCE_DIR}/cmake/dependencies.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/make_static_target.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/warnings.cmake)

include(CheckLanguage)
check_language(CUDA)
enable_language(CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_ARCHITECTURES 75)
set(CMAKE_CUDA_STANDARD ${CMAKE_CUDA_STANDARD_DEFAULT})
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# GLU library
# -----------
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "src/*.cpp" "src/*.c" "src/*.cu")
add_library(GLU ${sources})

# Alias for projects including GLU
add_library(GLU::GLU ALIAS GLU)

# Include directories
target_include_directories(GLU PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})


# Dependencies
target_link_libraries(GLU PRIVATE nicslu)

# Properties
set_target_properties(GLU
                      PROPERTIES POSITION_INDEPENDENT_CODE
                                 ON)
set_target_properties(GLU PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
make_target_static(GLU)
target_add_warnings(GLU)

# GLU command line interface
# --------------------------

add_executable(GLU-CLI "src/preprocess.c" "main.cpp")

# Dependencies
target_link_libraries(GLU-CLI PRIVATE GLU)
target_link_libraries(GLU-CLI PRIVATE nicslu)

# Properties
make_target_static(GLU-CLI)
target_add_warnings(GLU-CLI)
set_target_properties(GLU-CLI PROPERTIES OUTPUT_NAME "glu-cli")