cmake_minimum_required(VERSION 3.13)

include($ENV{DEVKIT_LOCATION}/cmake/pluginEntry.cmake)

set(PROJECT_NAME rpq9MultiSoftModDeformer)
project(${PROJECT_NAME} LANGUAGES C CXX)

set(MEL_FILE "${CMAKE_SOURCE_DIR}/src/AErpq9MultiSoftModTemplate.mel")
set(KERNEL_FILE "${CMAKE_SOURCE_DIR}/src/rpq9MultiSoftModDeformerKernels.cl")
set(GENERATED_DIR "${CMAKE_BINARY_DIR}/generated")
set(MEL_EMBEDDED_HDR "${GENERATED_DIR}/embeddedMel.h")
set(KERNEL_EMBEDDED_HDR "${GENERATED_DIR}/embeddedKernel.h")

add_custom_command(
  OUTPUT  "${MEL_EMBEDDED_HDR}"
  COMMAND "${CMAKE_COMMAND}" -E echo "[embedMel] generating: ${MEL_EMBEDDED_HDR}"
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${GENERATED_DIR}"
  COMMAND "${CMAKE_COMMAND}"
          -DMEL_FILE:PATH="${MEL_FILE}"
          -DOUT_FILE:PATH="${MEL_EMBEDDED_HDR}"
          -P "${CMAKE_SOURCE_DIR}/cmake/embedMel.cmake"
  DEPENDS "${MEL_FILE}" "${CMAKE_SOURCE_DIR}/cmake/embedMel.cmake"
)

add_custom_command(
  OUTPUT  "${KERNEL_EMBEDDED_HDR}"
  COMMAND "${CMAKE_COMMAND}" -E echo "[embedKernel] generating: ${KERNEL_EMBEDDED_HDR}"
  COMMAND "${CMAKE_COMMAND}" -E make_directory "${GENERATED_DIR}"
  COMMAND "${CMAKE_COMMAND}"
          -DKERNEL_FILE:PATH="${KERNEL_FILE}"
          -DOUT_FILE:PATH="${KERNEL_EMBEDDED_HDR}"
          -P "${CMAKE_SOURCE_DIR}/cmake/embedKernel.cmake"
  DEPENDS "${KERNEL_FILE}" "${CMAKE_SOURCE_DIR}/cmake/embedKernel.cmake"
)

set(RESOURCES_FILES 
    "${KERNEL_FILE}"
    "${MEL_FILE}"
    )

file(GLOB SRC_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB HED_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.h")

set(SOURCE_FILES
    ${SRC_FILES}
    ${HED_FILES}
   ${RESOURCES_FILES}
)

set(LIBRARIES
     OpenMaya
     OpenMayaAnim
     OpenMayaRender
     Foundation
     clew
)

if(MSVC)
    add_compile_options(/utf-8)
else()
    add_compile_options(-finput-charset=UTF-8 -fexec-charset=UTF-8)
endif()

find_tbb()

build_plugin()

# Put OUTPUT in the build graph
add_custom_target(embedAll DEPENDS 
                "${MEL_EMBEDDED_HDR}"
                "${KERNEL_EMBEDDED_HDR}")
add_dependencies(${PROJECT_NAME} embedAll)
set_source_files_properties("${MEL_EMBEDDED_HDR}"
                            "${KERNEL_EMBEDDED_HDR}"
                            PROPERTIES GENERATED TRUE)

target_include_directories(${PROJECT_NAME} PRIVATE "${GENERATED_DIR}")
target_sources(${PROJECT_NAME} PRIVATE 
                "${MEL_EMBEDDED_HDR}"
                "${KERNEL_EMBEDDED_HDR}")