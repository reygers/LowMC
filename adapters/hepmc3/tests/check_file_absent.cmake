if(NOT DEFINED CHECK_PATH)
    message(FATAL_ERROR "CHECK_PATH is required")
endif()

if(EXISTS "${CHECK_PATH}")
    message(FATAL_ERROR "Expected file to be absent: ${CHECK_PATH}")
endif()
