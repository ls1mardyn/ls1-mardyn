#!/bin/bash

while true
do
    ps axwu | grep di56giq5
    echo ""
    grep Huge /proc/meminfo
    echo ""
    sleep 3
done
