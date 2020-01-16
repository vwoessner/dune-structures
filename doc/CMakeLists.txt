add_document(
    TARGET readme
    OUTPUT_FILE readme.pdf
    SOURCES README.md
    PANDOC_DIRECTIVES --lua-filter ${CMAKE_CURRENT_SOURCE_DIR}/gitlab-math.lua
)
