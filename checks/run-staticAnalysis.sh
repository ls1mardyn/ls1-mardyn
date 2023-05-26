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
  if grep -q "using namespace std;" "$file"; then
    echo "\"using Log::global_log;\" was used in $file"
    warnings+="Do not use \"using Log::global_log;\"\n"
  fi

  # Check if file is UTF-8 (or ASCII)
  encoding=$(file -b --mime-encoding "$file")
  if [[ $encoding != "us-ascii" && $encoding != "utf-8" ]]; then
    # Fix
    #tmp_file=$(mktemp)
    #iconv -f "$encoding" -t utf-8 -o "$tmp_file" "$file"
    #mv -f "$tmp_file" "$file"
    echo "The following file is not ASCII/UTF-8 encoded: $file"
    warnings+="At least one file is not ASCII/UTF-8 encoded\n"
  fi

  # Check if file contains CRLF line endings
  fileinfo=$(file -k "$file")
  if [[ $fileinfo == *"CRLF"* ]]; then
    # Fix
    #sed -i 's/\r$//' "$file"
    echo "The following file contains CRLF line endings: $file"
    warnings+="At least one file contains CRLF line endings\n"
  fi

  # Check if file starts with BOM
  if [[ $fileinfo == *"BOM"* ]]; then
    # Fix
    #sed -i '1s/^\xEF\xBB\xBF//' "$file"
    echo "The following file starts with BOM: $file"
    warnings+="At least one file starts with BOM\n"
  fi

  # Check if file ends with newline
  if [[ -n "$(tail -c 1 "$file")" ]]; then
    # Fix
    #sed -i -e '$a\' "$file"
    echo "The following file does not end with newline: $file"
    warnings+="At least one file does not end with newline\n"
  fi

done

# Only print warnings once
warnings=$(printf "$warnings" | sort | uniq)
warnings="Warnings:\n"$warnings"\n\n"

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
    git fetch
    git switch master &> /dev/null
  else
    git switch ${currentVersion} &> /dev/null
  fi

  cpplint --filter=-whitespace/tab --linelength=200 --counting=detailed ${codeFiles} &> $rootFolder/staticAnalysis_${VERSION}.log

  grep "Category\|Total errors found" $rootFolder/staticAnalysis_${VERSION}.log > $rootFolder/staticAnalysis_${VERSION}_summary.log

done

git switch ${currentVersion} &> /dev/null

warnings+="\ncpplint:\n   New or fixed warnings/errors (master <-> new commit):\nmaster $(printf '%54s' " ") | new\n"

warnings+=$(diff -y $rootFolder/staticAnalysis_master_summary.log $rootFolder/staticAnalysis_new_summary.log)

tail $rootFolder/staticAnalysis*.log

printf "$warnings\n"
printf "$warnings\n" >> $GITHUB_STEP_SUMMARY
