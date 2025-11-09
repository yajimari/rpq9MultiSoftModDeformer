if(NOT DEFINED MEL_FILE)
  message(FATAL_ERROR "MEL_FILE not set")
endif()
if(NOT DEFINED OUT_FILE)
  message(FATAL_ERROR "OUT_FILE not set")
endif()

# Read MEL (UTF-8 without BOM)
file(READ "${MEL_FILE}" MEL_CODE_RAW)

# Write as header
file(WRITE "${OUT_FILE}" "/* Auto-generated from ${MEL_FILE} */\n#pragma once\n\n")
file(APPEND "${OUT_FILE}" "static const char* kEmbeddedMel = R\"EMB(\n")
file(APPEND "${OUT_FILE}" "${MEL_CODE_RAW}")
file(APPEND "${OUT_FILE}" "\n)EMB\";\n")