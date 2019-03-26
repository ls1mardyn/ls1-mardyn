#!/bin/bash
FILES=$(find ./src/ -name "*.cpp")
FILES+=" "
FILES+=$(find ./src/ -name "*.h")
FILES+=" "
FILES+=$(find ./src/ -name "*.hpp")
FILES+=" "
FILES+=$(find ./src/ -name "CMakeLists.txt")
echo "FILES=" $FILES

set -x
tar cfz Mardyn-src.tar.gz tools/standalone-generators/ libs/ makefile/ src/Makefile cmake CMakeLists.txt LICENSE $FILES
