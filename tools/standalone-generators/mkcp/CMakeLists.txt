cmake_minimum_required(VERSION 3.3)
message(STATUS "cmake building mkcp")

add_executable(mkcp main.cpp Domain.cpp Graphit.cpp Random.cpp)
target_include_directories(mkcp PRIVATE ${PROJECT_SOURCE_DIR})

include(GNUInstallDirs)
install(TARGETS mkcp RUNTIME DESTINATION /usr/bin/)

