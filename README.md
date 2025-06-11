# ChPTdiagram
Program for the specification and evaluation of diagrams in chiral perturbation theory (ChPT).

This is ultimately based on the FORM procedures used by Hans Bijnens, although it has developed to the point where no Hans-written code remains.
The laborious hand-definition of diagrams is replaced by code generation, and the complicated `pickout` procedure has been replaced by the FORM builtin `id,all`.

The Python program ChPTdiagram.py generates FORM code that can then be further customized using the procedures loaded under ChPTdiagram.hf.
Assuming the command `ChPTdiagram` executes ChPTdiagram.py with python 3, basic usage is
```
$ ChPTdiagram --generate-form-main <diagrams.xpt>
``` 
where `diagrams.xpt` contains definitions of the diagrams. Run 
```
$ ChPTdiagram --help
```
for further instructions, and inspect the generated code for guidance on how to use the FORM library.

Requires python 3, FORM (https://www.nikhef.nl/~form/) and my Python library for permutations (https://github.com/mssjo/permute/).
