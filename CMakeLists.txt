# CMakeList.txt : CMake project for GLU, include source and define
# project specific logic here.
#
cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

string(TIMESTAMP TODAY "%Y%m%d")
string(TIMESTAMP YEAR "%Y")

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project(GLU VERSION ${TODAY} LANGUAGES C CXX)

 string(TOLOWER ${PROJECT_NAME} lcname)

 option(GLU_CLI "Compile CLI for testing" OFF)

include(${PROJECT_SOURCE_DIR}/cmake/dependencies.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/make_static_target.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/warnings.cmake)

include(CheckLanguage)
check_language(CUDA)
enable_language(CUDA)

include(GNUInstallDirs)
# Add PIC option
if(NOT BUILD_SHARED_LIBS)
  option(WITH_PIC "Build library with PIC" ON)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_ARCHITECTURES 75)
set(CMAKE_CUDA_STANDARD ${CMAKE_CUDA_STANDARD_DEFAULT})
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# GLU library
# -----------
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "src/*.cpp" "src/*.c" "src/*.cu")
add_library(${lcname} ${sources})

# Alias for projects including GLU
add_library(${lcname}::${lcname} ALIAS ${lcname})

# Dependencies
target_link_libraries(${lcname} PRIVATE nicslu)

# GNU Install Directories
install(TARGETS ${lcname}
        EXPORT ${lcname}-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT ${lcname}-targets
        FILE CMake${PROJECT_NAME}Config.cmake
        NAMESPACE ${lcname}::
        DESTINATION share/CMake${PROJECT_NAME}c/)

# Properties
set_target_properties(${lcname} PROPERTIES OUTPUT_NAME ${lcname})
set_target_properties(${lcname} PROPERTIES POSITION_INDEPENDENT_CODE ${WITH_PIC})
set_target_properties(${lcname} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
make_target_static(${lcname})
target_add_warnings(${lcname})

# Include directories
target_include_directories(${lcname} PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# Install include directory
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/")

if(GLU_CLI)
    # GLU command line interface
    # --------------------------
    add_executable(${lcname}-CLI "src/preprocess.c" "main.cpp")

    # Dependencies
    target_link_libraries(${lcname}-CLI PRIVATE ${lcname})
    target_link_libraries(${lcname}-CLI PRIVATE nicslu)

    # Properties
    make_target_static(${lcname}-CLI)
    target_add_warnings(${lcname}-CLI)
    set_target_properties(${lcname}-CLI PROPERTIES OUTPUT_NAME "glu-cli")
endif()