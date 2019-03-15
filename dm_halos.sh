#!/bin/bash

DMfile = $(python3 main.py)
for fil in *.txt;
do
  if ["$fil" -ne "$DMfile.txt" -a "$fil" -ne "${DMfile}_nounits.txt"]
  then
    mv "$fil" "test - $fil"
  fi
done
