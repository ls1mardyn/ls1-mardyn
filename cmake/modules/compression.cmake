# compression wrapper
option(ENABLE_COMPRESSION "Use compression wrapper" ${ENABLE_COMPRESSION_WRAPPER})
if(ENABLE_COMPRESSION)
    message(STATUS "Using Compression wrapper.")
    option(ENABLE_LZ4 "Use LZ4 in compression wrapper" ON)
    include(lz4)
else()
    message(STATUS "Compression wrapper disabled.")
endif()
