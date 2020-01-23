add_document(
    TARGET readme
    OUTPUT_FILE readme.pdf
    SOURCES README.md
    PANDOC_DIRECTIVES --lua-filter ${CMAKE_CURRENT_SOURCE_DIR}/gitlab-math.lua
)

set_target_properties(readme PROPERTIES EXCLUDE_FROM_ALL 1)
