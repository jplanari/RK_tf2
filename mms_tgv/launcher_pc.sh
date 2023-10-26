#!/bin/bash

#timesteps=("1.0" "1e-1" "1e-2" "1e-3" "1e-4")
#schemes=("paramEuler" "heunRK2" "heunRK3" "stdRK4" "ps4p7q")

timesteps=("1e-4")
schemes=("ps4p7q")


for scheme in "${schemes[@]}"; do
  for dt in "${timesteps[@]}"; do
    foldername="$scheme-$dt"
    mkdir $foldername
    cp mms_launcher.cfg $foldername
    cd $foldername
    
    echo "Param:double:_TimeStep:$dt" >> mms_launcher.cfg
    echo "Param:str:RKmethod:$scheme" >> mms_launcher.cfg

    tf-sim mms_launcher.cfg

    cd ..
  done
done
