#!/bin/bash
# Author 
# Adriaan Larmuseau - 1ste master computer wetenschappen
# 
# Compiler
# tested in an zsh environmen
#
# Function
# measure edge computation performance

echo -n "" > ./results_edges.txt

make perf_edges > /dev/null 2>&1
if [[ $? -eq 0 ]];then
    it=3
    while [ $it -le 10 ]
    do
        echo "Running -- 2 ** $it" 
        echo "$it" | ./perf_edges >> results_edges.txt
        let it=$it+1
    done
else
    echo " failed to make test  -- please look into this !"
fi

echo -n "" > ./results_distance.txt

make shortest_distance > /dev/null 2>&1
if [[ $? -eq 0 ]];then
# test 3 different matrices
if [[ -e G_128_1.data ]];then
    echo "testing speed of distance for 0.001"
   it=0
   while [ $it -lt 10 ]
   do 
        ./shortest_distance --timesteps --blocksize 32 G_128_1.data >> results_distance.txt
        let it=$it+1
   done
   echo "-------------------------" >> results_distance.txt
else
    echo " missing file G_128_1.data !! generate using generate_edges "
fi

if [[ -e G_128_5.data ]];then
    echo "testing speed of distance for 0.005"
   it=0
   while [ $it -lt 10 ]
   do 
        ./shortest_distance --timesteps --blocksize 32 G_128_5.data >> results_distance.txt
        let it=$it+1
   done
   echo "-------------------------" >> results_distance.txt
else
    echo " missing file G_128_5.data !! generate using generate_edges "
fi

if [[ -e G_128_25.data ]];then
    echo "testing speed of distance for 0.025"
   it=0
   while [ $it -lt 10 ]
   do 
        ./shortest_distance --timesteps --blocksize 32 G_128_25.data >> results_distance.txt
        let it=$it+1
   done
   echo "-------------------------" >> results_distance.txt
else
    echo " missing file G_128_25.data !! generate using generate_edges "
fi

else
    
    echo" failed to make ./shortest_distance !!"
fi
