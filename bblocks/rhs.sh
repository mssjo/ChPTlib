#!/bin/bash

# Creates the file describing how to distribute $1 powers on $2 building blocks
# (the right-hand side of the power distribution expression)

ipart --form -uf $2 $1 > bblocks/rhs/$1on$2.hf
for i in $(seq 1 $2)
do
    sed -i -E -e "s/ipart\(([0-9]+),?/bb$i(\1,?x$i) * ipart\(/"  bblocks/rhs/$1on$2.hf
done
sed -i -E -e "s/ \* ipart\(\)//"  bblocks/rhs/$1on$2.hf
cat bblocks/rhs/$1on$2.hf
