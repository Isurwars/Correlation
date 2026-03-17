# Correlation - Liquid and Amorphous Solid Analysis Tool

# Readers as an OBJECT library to ensure static registration triggers
file(GLOB READER_SOURCES CONFIGURE_DEPENDS
    "*Reader.cpp"
)

add_library(readers_obj OBJECT
    ${READER_SOURCES}
    ReaderFactory.cpp
)

target_include_directories(readers_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
