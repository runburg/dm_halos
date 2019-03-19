#!/usr/bin/env bash

DMfile="df_nfw"
for fil in *.txt; do
  if [[ "$fil" != "$DMfile.txt" && "$fil" != "${DMfile}_nounits.txt" ]]; then
    mv "$fil" "${DMfile}_$fil"
  fi
done

mkdir $DMfile/
mv *.txt $DMfile/
