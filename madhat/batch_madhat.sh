#!/usr/bin/env bash

path="/Users/runburg/github/dm_halos/madhat/"
beta="0.9500"
mv j_*.txt models/;
mv channel_*.txt channels/;
rm ./Output/*.out;
for file in models/*.txt; do
  echo $file;
  mv $file ./;
  for channel in channels/*.txt; do
    echo $channel
    mv $channel ./;
    chan=${channel##*/};
    ./MADHAT j_*.txt $beta channel_*.txt;
    mv Output/*.out jackoutput/${chan%.txt}_${file##*/};
    mv channel_*.txt channels/;
  done
  mv j_*.txt models/;
done
