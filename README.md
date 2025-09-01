# ChPTlib
Python program and FORM library for the specification and evaluation of diagrams in chiral perturbation theory (ChPT).

This is ultimately based on the FORM procedures used by Hans Bijnens, although it has developed to the point where no Hans-written code remains.
The laborious hand-definition of diagrams is replaced by code generation, and the complicated `pickout` procedure has been replaced by the FORM builtin `id,all`.

The Python program `ChPT.py` generates FORM code that can then be further customized using the FORM library, which may be loaded with `#include ChPTlib.hf`.
Assuming the command `ChPT` executes `ChPTpy` with python 3, basic usage is
```
$ ChPT --generate-form-main <diagrams.chpt>
``` 
where `diagrams.chpt` contains definitions of the diagrams. Run 
```
$ ChPT --help
```
for further instructions, and inspect the generated code (`diagrams.frm` in the example above) for guidance on how to use the FORM library.

Requires python 3, FORM (https://www.nikhef.nl/~form/) and my Python library for permutations (https://github.com/mssjo/permute/).
