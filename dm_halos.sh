#!/usr/bin/env bash


DMfile="fe_GC_NFW"
for fil in *.txt; do
  if [[ "$fil" != "$DMfile.txt" && "$fil" != "${DMfile}_nounits.txt" ]]; then
    mv "$fil" "$DMfile_$fil"
  fi
done

mkdir $DMfile/
mv *.txt $DMfile/
