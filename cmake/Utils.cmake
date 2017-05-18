# ###############################################################################################
# Builds a file from a template only if the file has changed
# Saving a new copy of a template triggers recompilation of its dependencies. It should be done
# only if required to prevent excessive recompilation.
# Usage:
#   build_template("${INCLUDE_PATH}/header_file.h.in" "${INCLUDE_PATH}/header_file.h")
macro(build_template TEMPLATE_FILE OUTPUT_FILE)
    configure_file(${TEMPLATE_FILE} "${OUTPUT_FILE}.temp")

    if(EXISTS ${OUTPUT_FILE})
        file(SHA256 ${OUTPUT_FILE} _last_hash)
        file(SHA256 "${OUTPUT_FILE}.temp" _current_hash)
        if(NOT _last_hash EQUAL _current_hash)
            file(RENAME "${OUTPUT_FILE}.temp" ${OUTPUT_FILE})
        endif()
        unset(_last_hash)
        unset(_current_hash)
    else()
        file(RENAME "${OUTPUT_FILE}.temp" ${OUTPUT_FILE})
    endif()

    file(REMOVE "${OUTPUT_FILE}.temp")
endmacro()