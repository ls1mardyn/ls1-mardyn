cmake_minimum_required(VERSION 3.3)
message(STATUS "cmake building mktcts")

add_executable(mktcts main.cpp Domain.cpp Random.cpp)
target_include_directories(mktcts PRIVATE ${PROJECT_SOURCE_DIR})

include(GNUInstallDirs)
install(TARGETS mktcts RUNTIME DESTINATION /usr/bin/)

