add_executable(main "main.cpp")
target_link_libraries(main PRIVATE ioha ioh)

add_executable(fGA "fGA.cpp")
target_link_libraries(fGA PRIVATE ioha ioh)