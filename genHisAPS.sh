#!/bin/bash
width=$1
skip=$2
for temp in 298.0 330.0 360.0;
  do 
  echo $temp
  cd $temp
  for mu in m*;
  do 
    echo $mu
    cd $mu
    hisFile.py -r 't'$temp''$mu'g' -w $width -o cassandra -s $skip
    cd ../
  done; 
  cd ..
done;
