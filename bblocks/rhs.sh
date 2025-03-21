#!/bin/bash

for i in $(seq 1 $2)
do
    ipart --form -uf $2 $1 |
        sed -E "s/ipart\(([0-9]+),?/bb$i(\1,?x$i) * ipart\(/" |
        sed -E "s/ \* ipart()//" > bblocks/rhs/$1on$2.hf
done
