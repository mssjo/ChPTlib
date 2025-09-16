include names.mk

SHELL := /bin/bash

FORM = form -d QUIET
TFORM = tform -w4 -l -d QUIET

$(VERTEXFILE): ChPT_lagrexpand.frm ChPT_p$(ORDER)lagrangian.hf ChPT_bblocks.hf ChPT_bbexpand.frm
	$(info Generating $(VERTEXFILE))
	$(FORM) $(FORMOPTS) ChPT_lagrexpand.frm

$(BBLOCKFILE) : ChPT_bblocks.hf ChPT_bbexpand.frm
	mkdir -p $(BBLOCKDIR)/$(BBLOCKNAME)
	$(info Generating $(BBLOCKFILE))
	$(FORM) $(FORMOPTS) -d BBLOCK=$(BBLOCK) -d RANK=$(RANK) ChPT_bbexpand.frm

$(BBLOCKDIR)/rhs/%on$(NBBLOCK).hf : $(BBLOCKDIR)/rhs.sh
	mkdir -p $(BBLOCKDIR)/rhs
	$(BBLOCKDIR)/rhs.sh $* $(NBBLOCK)

$(CHPTDIR)/flavs/flav$(NM).hf : make_flav.py
	$(info Generating flavs/flav$(NM).hf)
	mkdir -p flavs
	python make_flav.py $(NM)
$(CHPTDIR)/flavs/flav$(NM)_SUN.hf : make_flav.py
	$(info Generating flavs/flav$(NM)_SUN.hf)
	mkdir -p flavs
	python make_flav.py $(NM) SUN

.PHONY: clean, vertex, bblock
vertex : $(VERTEXFILE)
bblock : $(BBLOCKFILE)
clean :
	rm -f vertices/*
	rm -rf bblocks/bb*
	rm -f bblocks/rhs/*
	rm -f partitions/*
