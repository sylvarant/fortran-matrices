#!/bin/bash
# Author 
# Adriaan Larmuseau - 1ste master computer wetenschappen
# 
# Compiler
# tested in an zsh environmen
#
# Function
# run our experiment for std multiplication

cat /proc/cpuinfo > info.cpu

max=2048
maxit=10
commands=() #( --3loops --matmul --dgemm )

echo -n "" > ./results_std.txt
echo -n "" > ./results_z.txt

make generate > /dev/null 2>&1

i=7
while [ $i -le 11 ]
do
    echo "$i" | ./generate; let i=$i+1
done

for comm in ${commands[@]}
do
	echo "Running :: $comm"
	it=128
	while [ $it -le $max ]
	do
		echo "$comm :: $it" >> results_std.txt
        j=0
        while [ $j -le $maxit ] 
        do
		    ./matrix_product A_$it.data B_$it.data --time $comm >> results_std.txt
            let j=$j+1
        done
		let it=$it*2
	done
done

it=128
while [ $it -le $max ]
do
	bl=1
	echo "Running --matrix_product -> $it"
	while [ $bl -le $it ]
	do
		echo "$it :: $bl ">> results_z.txt 
        j=0
        while [ $j -le $maxit ] 
        do
		    ./matrix_product A_$it.data B_$it.data --time --zmatrix --blocksize $bl >> results_z.txt
            let j=$j+1
        done
		let bl=$bl*2
	done
	let it=$it*2
done
		





