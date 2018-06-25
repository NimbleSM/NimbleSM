#!/bin/bash

# bash shell script to backup work
# note: to unzip and extract, use "tar -xzf <file name>"
# to view contents, use "tar tfz <file name>"
# to unzip (without extracting), use gunzip

OF=NimbleSM_$(date +%y)_$(date +%m)_$(date +%d)_$(date +%H)$(date +%M).tgz

FILES="AUTHORS *.md *.txt ./cmake ./examples ./scripts ./src ./test"

echo
echo Creating backup at $OF ...

tar -czf $OF $FILES

echo
echo listing contents of tar file ...
echo

tar tfz $OF

echo
