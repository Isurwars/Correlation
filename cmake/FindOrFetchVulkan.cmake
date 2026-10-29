# -----------------------------------------------------------
# Find or Fetch Vulkan Development Headers and Dynamic Loader
# -----------------------------------------------------------

option(CORRELATION_ENABLE_VULKAN "Enable hardware-accelerated Vulkan 3D rendering backend" ON)

if(CORRELATION_ENABLE_VULKAN AND BUILD_GUI AND NOT (BUILD_WASM OR EMSCRIPTEN))
  # 1. Attempt system discovery
  find_package(Vulkan QUIET)

  set(_VULKAN_HEADERS_FOUND FALSE)
  if(Vulkan_FOUND AND Vulkan_INCLUDE_DIRS)
    foreach(_INC_DIR ${Vulkan_INCLUDE_DIRS})
      if(EXISTS "${_INC_DIR}/vulkan/vulkan.h")
        set(_VULKAN_HEADERS_FOUND TRUE)
        break()
      endif()
    endforeach()
  endif()

  # 2. Fetch Khronos Vulkan-Headers if system headers are missing
  if(NOT _VULKAN_HEADERS_FOUND)
    message(STATUS "System Vulkan headers not found. Fetching Khronos Vulkan-Headers via FetchContent...")
    correlation_push_shared_libs(OFF)
    FetchContent_Declare(
      vulkan_headers
      GIT_REPOSITORY https://github.com/KhronosGroup/Vulkan-Headers.git
      GIT_TAG v1.4.304
      GIT_SHALLOW TRUE
      SYSTEM
      EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(vulkan_headers)
    correlation_pop_shared_libs()
  endif()

  # 3. Locate Vulkan dynamic loader
  find_library(Vulkan_LOADER_LIB
    NAMES vulkan vulkan-1 libvulkan libvulkan.so.1
    HINTS /usr/lib /usr/lib64 /usr/local/lib /opt/homebrew/lib
  )

  # 4. Synthesize imported target Vulkan::Vulkan if needed
  if(NOT TARGET Vulkan::Vulkan)
    add_library(Vulkan::Vulkan INTERFACE IMPORTED GLOBAL)
    if(TARGET Vulkan::Headers)
      target_link_libraries(Vulkan::Vulkan INTERFACE Vulkan::Headers)
    elseif(Vulkan_INCLUDE_DIRS)
      target_include_directories(Vulkan::Vulkan INTERFACE ${Vulkan_INCLUDE_DIRS})
    endif()

    if(Vulkan_LOADER_LIB)
      target_link_libraries(Vulkan::Vulkan INTERFACE ${Vulkan_LOADER_LIB})
    elseif(Vulkan_LIBRARIES)
      target_link_libraries(Vulkan::Vulkan INTERFACE ${Vulkan_LIBRARIES})
    endif()
  endif()

  if(Vulkan_LOADER_LIB OR Vulkan_LIBRARIES)
    set(CORRELATION_HAS_VULKAN TRUE CACHE INTERNAL "Vulkan backend active")
    add_compile_definitions(CORRELATION_USE_VULKAN=1)
    message(STATUS "Vulkan 3D rendering backend ENABLED (Loader: ${Vulkan_LOADER_LIB})")
  else()
    set(CORRELATION_HAS_VULKAN FALSE CACHE INTERNAL "Vulkan backend active")
    message(WARNING "Vulkan loader library not found; Vulkan hardware acceleration disabled.")
  endif()

  # 5. Locate glslc shader compiler
  find_program(GLSLC_EXECUTABLE NAMES glslc HINTS ENV VULKAN_SDK PATH_SUFFIXES bin /usr/bin /usr/local/bin)
  if(GLSLC_EXECUTABLE)
    message(STATUS "Found glslc shader compiler: ${GLSLC_EXECUTABLE}")
  else()
    message(STATUS "glslc not found; using bundled pre-compiled SPIR-V C++ fallback headers.")
  endif()

else()
  set(CORRELATION_HAS_VULKAN FALSE CACHE INTERNAL "Vulkan backend active")
endif()
