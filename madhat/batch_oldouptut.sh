#!/usr/bin/env bash

path="/Users/runburg/github/dm_halos/madhat/"
beta="0.9500"
mv j_*.txt models/;
mv channel_*.txt channels/;
rm ./Output/*.out;

mv DwarfSets/Set_1.dat ./
for channel in channels/*.txt; do
  echo $channel
  mv $channel ./;
  chan=${channel##*/};
  ./MADHAT Set_1.dat $beta channel_*.txt;
  mv Output/*.out jackoutput/refpaperoutput/${chan%.txt}_j_s.txt;
  mv channel_*.txt channels/;
done
mv Set_1.dat DwarfSets/
mv DwarfSets/Set_5.dat ./
for channel in channels/*.txt; do
  echo $channel
  mv $channel ./;
  chan=${channel##*/};
  ./MADHAT Set_5.dat $beta channel_*.txt;
  mv Output/*.out jackoutput/refpaperoutput/${chan%.txt}_j_som.txt;
  mv channel_*.txt channels/;
done
mv Set_5.dat DwarfSets/
