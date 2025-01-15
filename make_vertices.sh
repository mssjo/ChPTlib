#!/bin/bash

set -e

runform () {
    if [ -z "${4// }" ]
    then
        echo "O(p^$1) $2-meson $3-photon vertex"
    else
        echo "O(p^$1) $2-meson $3-photon vertex with options $4"
    fi
    form -l -d NP=$1 -d NM=$2 -d NV=$3 $4 ChPTdiagram_lagrangian.frm > /dev/null
}

rm vertices/*

for nf in "" "-d NF=3" "-d NF=2 -d NOCAYHAM"
do
    for par in "" "-d GENPAR" "-d KALPAR"
    do
        for ord in 2 4 6
        do
            runform $ord 0 2 "$nf $par"
            runform $ord 2 0 "$nf $par"
            runform $ord 2 1 "$nf $par"
            runform $ord 2 2 "$nf $par"
        done
        runform 2 4 0 "$nf $par"
        runform 4 4 0 "$nf $par"
        runform 2 4 1 "$nf $par"
        runform 2 4 2 "$nf $par"
        runform 2 6 0 "$nf $par"
        runform 4 6 0 "$nf $par"
        runform 2 8 0 "$nf $par"
    done
done
