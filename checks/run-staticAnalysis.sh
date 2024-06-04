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
exit_code="0"

codeFiles=$( find $rootFolder/src -name "*.h" -or -name "*.cpp" )

# Similar to check_format of MegaMol repo
for file in $codeFiles; do

  # Check if using namespace std is used
  if grep -q "using namespace std;" "$file"; then
    echo "\"using namespace std;\" was used in $file"
    warnings+="- Do not use \"using namespace std;\"\n"
    exit_code="1"
  fi

  # Check if using Log::global_log is used
  if grep -q "using Log::global_log;" "$file"; then
    echo "\"using Log::global_log;\" was used in $file"
    warnings+="- Do not use \"using Log::global_log;\"\n"
    exit_code="1"
  fi

  # Check if file is UTF-8 (or ASCII)
  encoding=$(file -b --mime-encoding "$file")
  if [[ $encoding != "us-ascii" && $encoding != "utf-8" ]]; then
    echo "The following file is not ASCII/UTF-8 encoded: $file ($encoding)"
    echo "  Fix with:"
    echo "    tmp_file=\$(mktemp)"
    echo "    iconv -f \"\$(file -b --mime-encoding \"$file\")\" -t utf-8 -o \"\$tmp_file\" \"\$file\""
    echo "    mv -f \"\$tmp_file\" \"\$file\""
    warnings+="- At least one file is not ASCII/UTF-8 encoded\n"
    exit_code="1"
  fi

  # Check if file contains CRLF line endings
  fileinfo=$(file -k "$file")
  if [[ $fileinfo == *"CRLF"* ]]; then
    echo "The following file contains CRLF line endings: $file"
    echo "  Fix with:"
    echo "    sed -i 's/\r$//' \"$file\""
              sed -i 's/\r$//' "$file"
    warnings+="- At least one file contains CRLF line endings\n"
    exit_code="1"
  fi

  # Check if file starts with BOM
  if [[ $fileinfo == *"BOM"* ]]; then
    echo "The following file starts with BOM: $file"
    echo "  Fix with:"
    echo "    sed -i '1s/^\xEF\xBB\xBF//' \"$file\""
    warnings+="- At least one file starts with BOM\n"
    exit_code="1"
  fi

  # Check if file ends with newline
  if [[ -n "$(tail -c 1 "$file")" ]]; then
    echo "The following file does not end with newline: $file"
    echo "  Fix with:"
    echo "    sed -i -e '\$a\' \"$file\""
    warnings+="- At least one file does not end with newline\n"
    exit_code="1"
  fi

done

if [ "$exit_code" = "1" ]; then
    printf -- "\n\n"  # Some space to make output clearer
    
    # Only print warnings once to job summary
    warnings=$(printf -- "$warnings" | sort | uniq)
    warnings="# Warnings\n"$warnings"\n\n"

    printf -- "\n$warnings\n"  # Print to job output
    printf -- "\n$warnings\n" >> $GITHUB_STEP_SUMMARY  # Print to job summary
    printf -- "\nSee job step for details\n" >> $GITHUB_STEP_SUMMARY
else
    printf -- "\nNo warnings\n"  # Print to job output
    printf -- "\nNo warnings :rocket:\n" >> $GITHUB_STEP_SUMMARY  # Print to job summary
fi

exit $exit_code
