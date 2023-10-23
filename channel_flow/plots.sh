#!/bin/bash

simnames=("cfl" "sat")
schemes=("heunRK2")
res=("180")
nprobes=5

cd probes

for simname in "${simnames[@]}"; do
  for scheme in "${schemes[@]}"; do
    for re in "${res[@]}"; do
      foldername="$re-$scheme-$simname"
      re_num=$(awk "BEGIN { printf \"%.6f\", $re }")
      probebase="cflow-${simname}_test-Re_tau${re_num}-${scheme}"

      if ! test -d $foldername; then
        mkdir $foldername
      fi

      tf-probe "$probebase.tfpb"
      cp *dat $foldername
      rm *dat

      cd $foldername
 
      if ! test -d "ux"; then
        mkdir "ux"
      fi
      if ! test -d "uy"; then
        mkdir "uy"
      fi     
      if ! test -d "uz"; then
        mkdir "uz"
      fi
      if ! test -d "p"; then
        mkdir "p"
      fi
      
      for (( i=1; i<=$nprobes; i++ ))
      do
        echo $i
        python3 ../../plot_probes.py "$probebase.p$i.dat" $foldername $i   
      done
      cd .. 
    done    
  done
done
