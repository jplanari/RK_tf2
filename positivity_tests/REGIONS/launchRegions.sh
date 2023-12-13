#!/bin/bash

if [ ! -d "plots" ]; then
  mkdir "plots"
fi

for((i=1;i<=5;i++)); do
  python3 calculateRegions.py $i
  echo "Regions for RK$i completed."
done

