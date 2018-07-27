#!/bin/bash
FILES=$(find ./src/ -name "*.cpp")
FILES+=" "
FILES+=$(find ./src/ -name "*.h")
FILES+=" "
FILES+=$(find ./src/ -name "*.hpp")
echo "FILES=" $FILES

echo "tar cfz Mardyn-src.tar.gz tools/animake/ tools/mkCP/ tools/mkTcTS/ tools/mkesfera/  makefile/ \$FILES"
tar cfz Mardyn-src.tar.gz tools/animake/ tools/mkCP/ tools/mktcts/ tools/mkesfera/  makefile/ $FILES
