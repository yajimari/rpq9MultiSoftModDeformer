if(NOT DEFINED KERNEL_FILE)
  message(FATAL_ERROR "KERNEL_FILE not set")
endif()
if(NOT DEFINED OUT_FILE)
  message(FATAL_ERROR "OUT_FILE not set")
endif()

# Read Kernel (UTF-8 without BOM)
file(READ "${KERNEL_FILE}" KERNEL_CODE_RAW)

# Normalize line endings to avoid CRLF/LF differences
string(REPLACE "\r\n" "\n" KERNEL_CODE_RAW "${KERNEL_CODE_RAW}")
string(REPLACE "\r"   "\n" KERNEL_CODE_RAW "${KERNEL_CODE_RAW}")

# Split into smaller chunks so MSVC does not see one huge string literal
# Keep this comfortably below the compiler limit.
set(CHUNK_SIZE 8000)

string(LENGTH "${KERNEL_CODE_RAW}" KERNEL_LEN)

file(WRITE "${OUT_FILE}" "/* Auto-generated from ${KERNEL_FILE} */\n")
file(APPEND "${OUT_FILE}" "#pragma once\n\n")
file(APPEND "${OUT_FILE}" "static const char* kEmbeddedKernel =\n")

set(POS 0)
while(POS LESS KERNEL_LEN)
  math(EXPR REMAIN "${KERNEL_LEN} - ${POS}")

  if(REMAIN GREATER CHUNK_SIZE)
    set(CUR_SIZE ${CHUNK_SIZE})
  else()
    set(CUR_SIZE ${REMAIN})
  endif()

  string(SUBSTRING "${KERNEL_CODE_RAW}" ${POS} ${CUR_SIZE} KERNEL_CHUNK)

  # Emit adjacent raw string literals; C/C++ concatenates them automatically
  file(APPEND "${OUT_FILE}" "R\"EMB(\n${KERNEL_CHUNK})EMB\"\n")

  math(EXPR POS "${POS} + ${CUR_SIZE}")
endwhile()

file(APPEND "${OUT_FILE}" ";\n")