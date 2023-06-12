#!/bin/bash
# Script to run static code analysis

#set strict pipefail option
#set -eo pipefail

# A path to the root folder of ls1 mardyn can be set via an argument; if not, the parent folder of "checks" is used
if [ $# -eq 1 ]
then
  rootFolder=$1
else
  rootFolder=$PWD/..
fi

echo "Running in $rootFolder"

warnings=""

codeFiles=$( find $rootFolder/src -name "*.h" -or -name "*.cpp" -printf "%p " )

# Similar to check_format of MegaMol repo
for file in $codeFiles; do

  # Check if using namespace std is used
  if grep -q "using namespace std;" "$file"; then
    echo "\"using namespace std;\" was used in $file"
    warnings+="Do not use \"using namespace std;\"\n"
  fi

  # Check if using Log::global_log is used
  if grep -q "using Log::global_log;" "$file"; then
    echo "\"using Log::global_log;\" was used in $file"
    warnings+="Do not use \"using Log::global_log;\"\n"
  fi

  # Check if file is UTF-8 (or ASCII)
  encoding=$(file -b --mime-encoding "$file")
  if [[ $encoding != "us-ascii" && $encoding != "utf-8" ]]; then
    echo "The following file is not ASCII/UTF-8 encoded: $file ($encoding)"
    echo "  Fix with:"
    echo "    tmp_file=\$(mktemp)"
    echo "    iconv -f \"\$(file -b --mime-encoding \"$file\")\" -t utf-8 -o \"\$tmp_file\" \"\$file\""
    echo "    mv -f \"\$tmp_file\" \"\$file\""
    warnings+="At least one file is not ASCII/UTF-8 encoded\n"
  fi

  # Check if file contains CRLF line endings
  fileinfo=$(file -k "$file")
  if [[ $fileinfo == *"CRLF"* ]]; then
    echo "The following file contains CRLF line endings: $file"
    echo "  Fix with:"
    echo "    sed -i 's/\r$//' \"$file\""
              sed -i 's/\r$//' "$file"
    warnings+="At least one file contains CRLF line endings\n"
  fi

  # Check if file starts with BOM
  if [[ $fileinfo == *"BOM"* ]]; then
    echo "The following file starts with BOM: $file"
    echo "  Fix with:"
    echo "    sed -i '1s/^\xEF\xBB\xBF//' \"$file\""
    warnings+="At least one file starts with BOM\n"
  fi

  # Check if file ends with newline
  if [[ -n "$(tail -c 1 "$file")" ]]; then
    echo "The following file does not end with newline: $file"
    echo "  Fix with:"
    echo "    sed -i -e '\$a\' \"$file\""
    warnings+="At least one file does not end with newline\n"
  fi

done

printf "\n\n\n"  # Some space to make output clearer

# Only print warnings once to job summary
warnings=$(printf "$warnings" | sort | uniq)
warnings="# Warnings\n"$warnings"\n\n"

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

currentVersion=$(git rev-parse --abbrev-ref HEAD)

for VERSION in "new" "master"
do
  echo "Checking ${VERSION} version"

  if [[ "${VERSION}" == "master" ]]
  then
    git fetch &> /dev/null
    git switch master &> /dev/null
  else
    git switch ${currentVersion} &> /dev/null
  fi

  cpplint --filter=-whitespace/tab --linelength=200 --counting=detailed ${codeFiles} &> $rootFolder/staticAnalysis_${VERSION}.log

  grep "Category\|Total errors found" $rootFolder/staticAnalysis_${VERSION}.log > $rootFolder/staticAnalysis_${VERSION}_summary.log

done

git switch ${currentVersion} &> /dev/null

# Use code block for monospace font in Markdown
warnings+="\n# cpplint\n"
warnings+=" For details run \`cpplint --filter=-whitespace/tab --linelength=200 --counting=detailed \$( find src -name \"*.h\" -or -name \"*.cpp\" -printf \"%%p \" )\`"
warnings+=" New or fixed warnings/errors (master <-> new commit):\n"
warnings+="\`\`\`\n"
warnings+="master $(printf '%54s' " ") | new\n"

# Delete "--suppress-common-lines" to see all errors/warnings
warnings_cpplint=$(diff -y --suppress-common-lines $rootFolder/staticAnalysis_master_summary.log $rootFolder/staticAnalysis_new_summary.log)
warnings+=$warnings_cpplint

# Counts the categories in which new errors were introduced
exitcode=$(printf "$warnings_cpplint" | grep "Category" | awk '{if ($11 > $5) {count++}} END {print count}')

warnings+="\n\`\`\`\n"  # Close code block for monospace font
printf "\n$warnings\n"  # Print to job output
printf "\n$warnings\n" >> $GITHUB_STEP_SUMMARY  # Print to job summary

exit $exitcode
