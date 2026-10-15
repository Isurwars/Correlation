# SPDX-License-Identifier: AGPL-3.0-only
# Strict warning policy applied per-target (never globally, to keep third-party code quiet).

option(CORRELATION_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

add_library(correlation_warnings INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(correlation_warnings INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Wshadow>
        $<$<AND:$<BOOL:${CORRELATION_WARNINGS_AS_ERRORS}>,$<COMPILE_LANGUAGE:CXX>>:-Werror>)
elseif(MSVC)
    target_compile_options(correlation_warnings INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:/W4>
        $<$<AND:$<BOOL:${CORRELATION_WARNINGS_AS_ERRORS}>,$<COMPILE_LANGUAGE:CXX>>:/WX>)
endif()

# Link the warning policy into every first-party compile unit target.
function(correlation_enable_warnings)
    foreach(tgt IN LISTS ARGN)
        if(TARGET ${tgt})
            target_link_libraries(${tgt} PRIVATE correlation_warnings)
        endif()
    endforeach()
endfunction()
