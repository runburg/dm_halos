#!/bin/bash

DMfile = "fe_GC_NFW"
for fil in *.txt;
do
  if ["$fil" -ne "$DMfile.txt" -a "$fil" -ne "${DMfile}_nounits.txt"]
  then
    mv "$fil" "test - $fil"
  fi
done

mkdir $DMfile/
mv *.txt $DMfile/
