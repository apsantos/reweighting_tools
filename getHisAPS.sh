#!/bin/bash
run=$1
make=$2
for TEMP in 298.0 330.0 360.0;
do
 cd ../$TEMP/
  for mu in m*;
  do 
    cd $mu
    cp 'hist'$TEMP''$mu'g'$run'.dat' ../../reweigh/
    cd ../
  done
done
cd ../reweigh
if [ "$make" -eq "1" ]
  then
  echo '1' > 'input_hs.dat' 
  echo $run >> 'input_hs.dat' 
  for his in his*.dat;
  do 
    NAME=`echo "$his" | cut -d"s" -f2`
    NAME=`echo "$NAME" | cut -d"$run" -f1`
    echo $NAME >> 'input_hs.dat'
  done
fi
