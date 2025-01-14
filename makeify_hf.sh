#!/bin/bash

# Usage: ./makeify_hf.sh FILE.hf 
# Outputs FILE.make

# Script that converts a FORM header file into an equivalent Make header file,
# assuming it contains only simple definitions and ifdefs,
# using the fact that that subset of their syntaxes are rather similar.

# Used here to generate names.make from names.hf so that filenames and tags stay consistent.

bt='`'

{
    echo "# Generated from $1.hf using makeify_hf.sh";
    echo "# DO NOT EDIT";
    echo "";
    cat $1.hf |
        sed -E \
            -e "s/$bt([A-Z\$()]+)'/\$(\1)/g" \
            -e 's/#(ifn?def) \$\(([A-Z]+)\)/\1 \2/' \
            -e 's/#(else|endif)/\1/' \
            -e 's/#include-? ([A-Za-z_]+).hf/include \1.make/' \
            -e 's/#define ([A-Z]+) "([^"]*)"/\1 = \2/' \
            -e 's/#redefine ([A-Z]+) "([^"]*)"/\1 := \2/' \
            -e 's/^\*/#/'
} | tee $1.make
