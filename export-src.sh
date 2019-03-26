#!/bin/bash
FILES=$(find ./src/ -name "*.cpp")
FILES+=" "
FILES+=$(find ./src/ -name "*.h")
FILES+=" "
FILES+=$(find ./src/ -name "*.hpp")
echo "FILES=" $FILES

set -x
tar cfz Mardyn-src.tar.gz tools/standalone-generators/ libs/armadillo/ libs/rapidxml/ makefile/ src/Makefile $FILES
