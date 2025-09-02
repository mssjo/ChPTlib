#!/bin/bash

printf "#define LOCALDIR \".\"\n#define CHPTDIR \"$(pwd)\"\n" > dirs.hf
printf "#!/bin/bash\npython3 $(pwd)/ChPT.py \$@\n" > ~/bin/ChPT
chmod +x ~/bin/ChPT
./makeify_hf.sh dirs > /dev/null
./makeify_hf.sh names > /dev/null

printf "NOTE: to complete the installation, add $(pwd) to FORMPATH\n"
