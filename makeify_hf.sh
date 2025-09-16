#!/bin/bash

bt='`'

{
    echo "# Generated from $1.hf using makeify_hf.sh";
    echo "# DO NOT EDIT";
    echo "";
    cat $1.hf |
        sed -E \
            -e ":x; s/$bt([A-Z\\\$()]+)'/\$(\1)/g; tx" \
            -e 's/#if *(\$\([A-Z]+\)) *== *([-+]?[0-9]+)/ifeq (\1,\2)/' \
            -e 's/#if *(\$\([A-Z]+\)) *!= *([-+]?[0-9]+)/ifneq (\1,\2)/' \
            -e 's/#if *(\$\([A-Z]+\))/ifneq (\1,0)/' \
            -e 's/#(ifn?def) \$\(([A-Z]+)\)/\1 \2/' \
            -e 's/#(else|endif)/\1/' \
            -e 's/#include-? ([A-Za-z_]+).hf/include \1.mk/' \
            -e 's/#define ([A-Z]+) "([^"]*)"/\1 = \2/' \
            -e 's/#redefine ([A-Z]+) "([^"]*)"/\1 := \2/' \
            -e 's/#define ([A-Z]+)/\1 = /' \
            -e 's/^\*/#/'
} | tee $1.mk
