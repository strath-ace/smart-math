cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# Change this section to correspond with your project info.
project (
  smart-math
  VERSION 0.1.0
  LANGUAGES CXX
)

set(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem ")

set(PROJECT_SUMMARY                            "Strathclyde Mechanical and Aerospace Researt Toolbox for Mathematics (SMART-MATH)")
set(${PROJECT_NAME}_VENDOR_NAME                "spaceART")
set(${PROJECT_NAME}_VENDOR_CONTACT             "spaceART@strathclyde.ac.uk")

set(PROJECT_PATH                               "${PROJECT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}     "${PROJECT_PATH}/cmake/Modules")
set(INCLUDE_PATH                               "${PROJECT_PATH}/include")
set(EIGEN_INCLUDE_PATH                         "${PROJECT_PATH}/include/LinearAlgebra")
set(SRC_PATH                                   "${PROJECT_PATH}/src")
set(TEST_SRC_PATH                              "${PROJECT_PATH}/test")
if(NOT EXTERNAL_PATH)
  set(EXTERNAL_PATH                            "${PROJECT_PATH}/external")
endif()
if(NOT DOCS_PATH)
  set(DOCS_PATH                                "${PROJECT_PATH}/docs")
endif()
set(LIB_PATH                                   "${PROJECT_BINARY_DIR}/lib")
set(LIB_NAME                                   "smart-math")
set(BIN_PATH                                   "${PROJECT_BINARY_DIR}/bin")
set(MAIN_NAME                                  "main")
set(TEST_PATH                                  "${PROJECT_BINARY_DIR}/test")
set(TEST_NAME                                  "test")
set(EXAMPLES_PATH                              "${PROJECT_BINARY_DIR}/examples")
set(EXAMPLES_NAME                              "examples")

option(BUILD_STATIC			       "Build static library"		OFF)
option(BUILD_MAIN                              "Build main function"            ON)
option(BUILD_DOCS                              "Build docs"                     OFF)
option(BUILD_TESTS                             "Build tests"                    OFF)
option(BUILD_EXAMPLES                          "Build examples"                 OFF)

include("cmake/Utils.cmake")
build_template("${INCLUDE_PATH}/smartmath.h.in" "${INCLUDE_PATH}/smartmath.h")

file(GLOB SRC
  "${SRC_PATH}/Utils/*.cpp"
)

set(MAIN_SRC
  "${SRC_PATH}/main.cpp"
)

file(GLOB TEST_SRC
  "${TEST_SRC_PATH}/LinearAlgebra/*.cpp"
  "${TEST_SRC_PATH}/Dynamics/*.cpp"
  "${TEST_SRC_PATH}/Events/*.cpp"
  "${TEST_SRC_PATH}/Integrators/*.cpp"
  "${TEST_SRC_PATH}/Utils/*.cpp"
)

file(GLOB_RECURSE LibFiles
    "${INCLUDE_PATH}/*.h"
    )
add_custom_target(headers SOURCES ${LibFiles})

# Set CMake build-type. If it not supplied by the user, the default built type is "Release".
if((NOT CMAKE_BUILD_TYPE)
  OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    message(STATUS "WARNING: building release version!")
    set(CMAKE_BUILD_TYPE Release)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "WARNING: building debug version!")
endif()

# Set platform-specific compiler flags.
if(WIN32)
  if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS -std=c++11")
    # We need to force MSVC to use static runtime:
    # (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace).
    foreach(flag_var
      CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      # Find all dynamic runtime (MD) references and replace with static (MT)
      if(${flag_var} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
  else(MSVC)
    set(CMAKE_LEGACY_CYGWIN_WIN32 0)
    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -Woverloaded-virtual -Wnon-virtual-dtor -std=c++11")
  endif()
else(WIN32)
  set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -Woverloaded-virtual -Wnon-virtual-dtor -std=c++11")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3 -fopenmp -std=c++11")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -std=c++11")
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage -std=c++11")
endif()

include_directories(${EIGEN_INCLUDE_PATH})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_PATH})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_PATH})
if(BUILD_STATIC) 
   add_library(${LIB_NAME} STATIC ${SRC})
else()
   add_library(${LIB_NAME} SHARED ${SRC})
endif()

set_target_properties(${LIB_NAME} PROPERTIES
	  VERSION ${PROJECT_VERSION}
	  SOVERSION ${PROJECT_VERSION}
)

if(BUILD_MAIN)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_PATH})
  add_executable(${MAIN_NAME} ${MAIN_SRC})
  target_link_libraries(${MAIN_NAME} ${LIB_NAME} ${MANDATORY_LIBRARIES})
endif()

if(BUILD_DOCS)
  find_package(Doxygen)

  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR
            "Doxygen is needed to build the documentation. Please make sure it is installed.")
  endif(NOT DOXYGEN_FOUND)

  configure_file(${PROJECT_PATH}/Doxyfile.in
                 ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

  # Add a custom target to run Doxygen whenever the project is built.
  # If you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the command below.
  add_custom_target(docs ALL COMMAND ${DOXYGEN_EXECUTABLE}
                    ${PROJECT_BINARY_DIR}/Doxyfile
                    SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
endif()

if(BUILD_TESTS)
  enable_testing()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_PATH})

  add_executable(${TEST_NAME} ${TEST_SRC})
  target_link_libraries(${TEST_NAME} ${LIB_NAME} ${MANDATORY_LIBRARIES})
  add_test(NAME ${TEST_NAME} COMMAND "${TEST_PATH}/${TEST_NAME}")

endif()

if(BUILD_EXAMPLES)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXAMPLES_PATH})
	add_subdirectory("${PROJECT_PATH}/examples")
endif()

# Install header files and library.
# Destination is set by CMAKE_INSTALL_PREFIX and defaults to usual locations, unless overridden by
# user.
install(DIRECTORY ${INCLUDE_PATH}/
        DESTINATION ${PROJECT_NAME}/include)
	#FILES_MATCHING PATTERN "*.h")
install(TARGETS ${LIB_NAME} DESTINATION ${PROJECT_NAME}/lib)
install(TARGETS ${BIN_NAME} DESTINATION ${PROJECT_NAME}/bin)

# Set up packager.
if(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
  include(InstallRequiredSystemLibraries)
  
  set(CPACK_SET_DESTDIR "on")
  set(CPACK_PACKAGING_INSTALL_PREFIX "/tmp")
  
  set(CPACK_PACKAGE_DESCRIPTION "${PROJECT_NAME}")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_NAME}")
  set(CPACK_PACKAGE_VENDOR "${PROJECT_NAME}_VENDOR_NAME")
  set(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
  set(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
  set(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
  set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_PATH}/README.md")
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_PATH}/LICENSE.md")
  set(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
  set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
  set(CPACK_GENERATOR "ZIP")

  if (UNIX)
    set(CPACK_GENERATOR "${CPACK_GENERATOR};DEB")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "strathclydesmart.slack.com")
    set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
    set(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})
  elseif (WIN32)
    # There is a bug in NSI that does not handle full unix paths properly. Make
    # sure there is at least one set of four (4) backlasshes.
    set(CPACK_GENERATOR "${CPACK_GENERATOR};NSIS")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
    set(CPACK_NSIS_CONTACT "${PROJECT_NAME}_VENDOR_CONTACT")
    set(CPACK_NSIS_MODIFY_PATH ON)
  endif()

  include(CPack)
endif()
