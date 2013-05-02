cmake_minimum_required (VERSION 2.6)
project (rst_all)

link_directories(/usr/lib)
add_subdirectory(librst)
add_subdirectory(projects/EmptyTab)
add_subdirectory(projects/TemplateTab)
