# SPDX-License-Identifier: AGPL-3.0-only
# Strict warning policy applied per-target (never globally, to keep third-party code quiet).

option(CORRELATION_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

add_library(correlation_warnings INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(correlation_warnings INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Wshadow -Wfloat-conversion -Wimplicit-float-conversion -Wshorten-64-to-32 -Werror=float-conversion -Werror=implicit-float-conversion -Werror=shorten-64-to-32>
        $<$<AND:$<BOOL:${CORRELATION_WARNINGS_AS_ERRORS}>,$<COMPILE_LANGUAGE:CXX>>:-Werror>)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(correlation_warnings INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Wshadow -Wfloat-conversion -Werror=float-conversion>
        $<$<AND:$<BOOL:${CORRELATION_WARNINGS_AS_ERRORS}>,$<COMPILE_LANGUAGE:CXX>>:-Werror>)
elseif(MSVC)
    target_compile_options(correlation_warnings INTERFACE
        $<$<COMPILE_LANGUAGE:CXX>:/W4 /wd4324 /we4244 /we4267 /we4305>
        $<$<AND:$<BOOL:${CORRELATION_WARNINGS_AS_ERRORS}>,$<COMPILE_LANGUAGE:CXX>>:/WX>)
    target_compile_definitions(correlation_warnings INTERFACE
        _CRT_SECURE_NO_WARNINGS)
endif()

# Link the warning policy into every first-party compile unit target.
function(correlation_enable_warnings)
    foreach(tgt IN LISTS ARGN)
        if(TARGET ${tgt})
            target_link_libraries(${tgt} PRIVATE correlation_warnings)
        endif()
    endforeach()
endfunction()
