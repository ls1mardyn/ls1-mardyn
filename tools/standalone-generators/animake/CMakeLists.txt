cmake_minimum_required(VERSION 3.3)
message(STATUS "cmake building animake")

add_executable(animake main.cpp Domain.cpp Random.cpp)
target_include_directories(animake PRIVATE ${PROJECT_SOURCE_DIR})

include(GNUInstallDirs)
install(TARGETS animake RUNTIME DESTINATION /usr/bin/)

