include names.make

SHELL := /bin/bash

FORM = tform -w4 -l -d QUIET
ifdef NP
	FORM := $(FORM) -d ORDER=$(NP)
endif


$(VERTEXFILE): ChPTdiagram_lagrexpand.frm ChPTdiagram_p$(NP)lagrangian.hf ChPTdiagram_bblocks.hf ChPTdiagram_bbexpand.frm
	$(info Generating $(VERTEXFILE))
	$(FORM) $(FORMOPTS) ChPTdiagram_lagrexpand.frm
# ifdef GENPAR
# 	sed -i -E ':x; s/a\(([0-9]+),?([0-9,]*)\)/a(\2)\1/g; tx; s/a\(\)/a/' $(VERTEXFILE)
# endif

$(BBLOCKFILE) : ChPTdiagram_bblocks.hf ChPTdiagram_bbexpand.frm
	mkdir -p $(BBLOCKDIR)/$(BBLOCKNAME)
	$(info Generating $(BBLOCKFILE))
	$(FORM) $(VERTEXOPTS) -d BBLOCK=$(BBLOCK) -d RANK=$(RANK) ChPTdiagram_bbexpand.frm

$(BBLOCKDIR)/rhs/%on$(NB).hf : partitions/%uf$(NB).hf $(BBLOCKDIR)/rhs.sh
	mkdir -p $(BBLOCKDIR)/rhs
	$(BBLOCKDIR)/rhs.sh $* $(NB)

$(CHPTDIR)/flavs/flav$(NM).hf : make_flav.py
	mkdir -p flavs
	python make_flav.py $(NM)

.PHONY: clean, vertex, bblock
vertex : $(VERTEXFILE)
bblock : $(BBLOCKFILE)
clean :
	rm -f vertices/*
	rm -rf bblocks/bb*
	rm -f bblocks/rhs/*
	rm -f partitions/*
