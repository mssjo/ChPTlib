#!/bin/bash

file=bblocks/rhs/$1on$2.hf
cp ipart/$1f$2.hf $file

for i in $(seq 1 $2)
do
    sed -i -E "s/ipart\(([0-9]+),?/bb$i(\1,?x$i) * ipart\(/" $file
done
sed -i "s/ \* ipart()//" $file
