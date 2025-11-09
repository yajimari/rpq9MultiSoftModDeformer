if(NOT DEFINED KERNEL_FILE)
  message(FATAL_ERROR "KERNEL_FILE not set")
endif()
if(NOT DEFINED OUT_FILE)
  message(FATAL_ERROR "OUT_FILE not set")
endif()

# Read Kernel (UTF-8 without BOM)
file(READ "${KERNEL_FILE}" KERNEL_CODE_RAW)

# Write as header
file(WRITE "${OUT_FILE}" "/* Auto-generated from ${KERNEL_FILE} */\n#pragma once\n\n")
file(APPEND "${OUT_FILE}" "static const char* kEmbeddedKernel = R\"EMB(\n")
file(APPEND "${OUT_FILE}" "${KERNEL_CODE_RAW}")
file(APPEND "${OUT_FILE}" "\n)EMB\";\n")