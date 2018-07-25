#!/bin/bash
ls -d -1 $PWD/*sec1* > flist1.txt
ls -d -1 $PWD/*sec2* > flist2.txt

python3 ../../../../../tools/ppls1/scripts/ppcoal_01.py

