cmake_minimum_required(VERSION 3.3)
message(STATUS "cmake building mkesfera")

add_executable(mkesfera main.cpp Domain.cpp Random.cpp)
target_include_directories(mkesfera PRIVATE ${PROJECT_SOURCE_DIR})

include(GNUInstallDirs)
install(TARGETS mkesfera RUNTIME DESTINATION /usr/bin/)

