#!/bin/bash
# Author 
# Adriaan Larmuseau - 1ste master computer wetenschappen
# 
# Compiler
# tested in an zsh environment, the echo to generate_edges is not
# portable so it is best to manually echo "9\n10" & "5\n 10" to
# ./generate_edges
#
# Function
# run all diagnostics

tests=( testconv testmult testread testmultedges test_distance test_steps )
programs=( matrix_product shortest_distance )
versionn=" Version 1.0 -- by Adriaan Larmuseau"

echo "== Running test Suite's , the quitter the better ! ==="

for t in ${tests[@]}
do
    make $t > /dev/null 2>&1
    if [[ $? -eq 0 ]];then
        echo "### Running test :: $t ###"
        ./$t
    else
        echo " failed to make test $t -- please look into this !"
    fi
done

# test additional paramters 
for p in ${programs[@]}
do
    echo "### Testing $p ###"
    make $p > /dev/null 2>&1 
    if [[ $? -eq 0 ]];then

        # should return nada
        ./$p 

        version=$(./$p --version)

        if [[ "$version" != "$versionn" ]]; then
            echo "Version command not working, for $p"
        fi

        ./$p --help > /dev/null 2>&1
        if [[ $? -lt 0 ]];then
            echo "help command failed for $p"
        fi

    else
        echo "$p is not compiling !"
    fi
done

# test some things unique to the different programs
make generate > /dev/null 2>&1
echo "### Testing basic command usage ###"
echo "9" | ./generate > /dev/null 2&>1
#echo "9\n10" | ./generate_edges
resulta=$(./matrix_product  A_512.data B_512.data --time) 
resultb=$(./shortest_distance G_512_10.data --time)
if [ $resulta = "0.0000000" ];then #|| 
    echo "-- !! time results not ok for matrix_product "
fi
if [ $resultb = "0.0000000" ]; then
    echo "-- !! time results not ok for shortest_distance "
fi
resultc=$(./matrix_product  A_512.data B_512.data --time --blocksize 8 ) 
resultd=$(./shortest_distance G_512_10.data --time --blocksize 8 )
#if [ "$resultc" -lt "$resulta" ];then 
#    echo "-- !! something wrong with blocksize in matrix_product"
#fi
# float's are not universal in shell just print
echo "$resultc $resultd $resulta $resultb"
echo "### Testing maxsteps ###"
./shortest_distance G_32_10.data --maxsteps 4 --steps | grep 5 
