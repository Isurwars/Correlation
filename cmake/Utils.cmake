# -----------------------------------------------------------
# Helper Function for Bundling DLLs on Windows
# -----------------------------------------------------------
function(bundle_windows_dlls target_name)
    if(WIN32)
        add_custom_command(
            TARGET ${target_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:${target_name}>
            $<TARGET_RUNTIME_DLLS:${target_name}>
            $<TARGET_FILE_DIR:${target_name}>
            COMMAND_EXPAND_LISTS
            COMMENT "Copying runtime DLLs for ${target_name}"
        )
    endif()
endfunction()

# -----------------------------------------------------------
# Helper: Link platform-specific dependencies for static Slint
# -----------------------------------------------------------
# When Slint is built from source as a static library, its Rust
# internals depend on system frameworks/libraries that are NOT
# propagated transitively through CMake.  Call this function on
# every target that links Slint::Slint.
function(link_slint_platform_deps target_name)
    if (APPLE)
        foreach(_fw IN ITEMS
            OpenGL CoreVideo Cocoa Carbon IOKit QuartzCore
            AppKit CoreGraphics Metal CoreFoundation Foundation Security)
            find_library(_${_fw}_FW ${_fw} REQUIRED)
            target_link_libraries(${target_name} PRIVATE ${_${_fw}_FW})
            unset(_${_fw}_FW CACHE)          # don't pollute the cache
        endforeach()
    elseif (NOT WIN32)
        find_package(Fontconfig REQUIRED)
        target_link_libraries(${target_name} PRIVATE ${Fontconfig_LIBRARIES})
    endif()
endfunction()

# Ensure that runtime DLLs are correctly included in the CPack installation
function(install_windows_dlls target_name dest_dir)
    if(WIN32)
        # Using a custom install script or directly installing the known TARGET_RUNTIME_DLLS
        # The generator expression $<TARGET_RUNTIME_DLLS:tgt> is only valid in CMake 3.21+
        install(FILES $<TARGET_RUNTIME_DLLS:${target_name}>
                DESTINATION ${dest_dir}
                COMPONENT correlation
                OPTIONAL)
    endif()
endfunction()
