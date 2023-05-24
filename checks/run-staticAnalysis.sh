#!/bin/bash
# Script to run static code analysis

codeFiles=$( find $PWD/../src -name "*.h" -or -name "*.cpp" )

# Check if using namespace std is introduced
if grep -q "using namespace std;" ${codeFiles}; then
	echo "Do not use \"using namespace std;\""
	return 1
fi


# Run cpplint

# Lines starting with one minus sign are filters which exclude certain rules

#-runtime/int,\
#-readability/alt_tokens,\
#-whitespace/blank_line,\
#-whitespace/braces,\
#-whitespace/comma,\
#-whitespace/comments,\
#-whitespace/indent,\
#-whitespace/operators,\
#-whitespace/parens,\

currentVersion=$(git rev-parse --verify HEAD)

for VERSION in "master" "new"
do
	if [[ "${VERSION}" == "master" ]]
	then
		git checkout origin/master
	else
		git checkout -
	fi

	cpplint --filter=\
	-whitespace/tab \
	--linelength=200 \
	--counting=detailed \
	${codeFiles} &> $PWD/staticAnalysis_${VERSION}.log

	grep "Category\|Total errors found" $PWD/staticAnalysis_${VERSION}.log | sort -k5 -nr > $PWD/staticAnalysis_${VERSION}_summary.log

done

echo "New or fixed warnings/errors (master <-> new commit):"
diff $PWD/staticAnalysis_master_summary.log $PWD/staticAnalysis_new_summary.log

