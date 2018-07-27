#!/bin/bash

find . -name "*.h" -o -name "*.c" -o -name "*.cpp" | xargs svn blame  | awk '{print $2}' | sort | uniq -c | sort -n -r
