# Generated from names.hf using makeify_hf.sh
# DO NOT EDIT

# Reads program options and defines preprocessor variables for
# - naming expressions, especially vertices
# - naming save files
# - calling external programs with the same options

include dirs.make

NAMESUFFIX = 

FORMOPTS = 
MAKEOPTS = 
ifdef NM
    ifdef NV
        ifdef NP
            HASINFO = 1
            FORMOPTS := -d NM=$(NM) -d NV=$(NV) -d NP=$(NP)
            MAKEOPTS := NM=$(NM) NV=$(NV) NP=$(NP)
        endif
    endif
endif

ifdef NFGENERAL
    MAKEOPTS := $(MAKEOPTS) NFGENERAL=1
else
    NAMESUFFIX := $(NAMESUFFIX)NF$(NF)
    FORMOPTS := $(FORMOPTS) -d NF=$(NF)
    MAKEOPTS := $(MAKEOPTS) NF=$(NF)
endif
ifdef CAYHAM
    NAMESUFFIX := $(NAMESUFFIX)_CAYHAM
    FORMOPTS := $(FORMOPTS) -d CAYHAM
    MAKEOPTS := $(MAKEOPTS) CAYHAM=1
endif
ifdef PAR
    NAMESUFFIX := $(NAMESUFFIX)_$(PAR)PAR
    FORMOPTS := $(FORMOPTS) -d PAR=$(PAR)
    MAKEOPTS := $(MAKEOPTS) PAR=$(PAR)
else
    NAMESUFFIX := $(NAMESUFFIX)_EXPPAR
    MAKEOPTS := $(MAKEOPTS) PAR=EXP
endif
ifdef TRF
    NAMESUFFIX := $(NAMESUFFIX)_TRF
    FORMOPTS := $(FORMOPTS) -d TRF
    MAKEOPTS := $(MAKEOPTS) TRF=1
endif
ifdef PKEMASS
    NAMESUFFIX := $(NAMESUFFIX)_PKE
    FORMOPTS := $(FORMOPTS) -d PKEMASS
    MAKEOPTS := $(MAKEOPTS) PKEMASS=1
endif
ifdef TRANSFORM
    NAMESUFFIX := $(NAMESUFFIX)_$(TRANSFORM)TRANSF
    FORMOPTS := $(FORMOPTS) -d TRANSFORM=$(TRANSFORM)
    MAKEOPTS := $(MAKEOPTS) TRANSFORM=$(TRANSFORM)
endif

MAKECMD = make --no-print-directory $(MAKEOPTS)
ifdef FORCEMAKE
    MAKECMD := $(MAKECMD) -B
endif

ifdef NAME
    SAVEDIR = $(LOCALDIR)/save
    SAVENAME = $(SAVEDIR)/$(NAME)$(NAMESUFFIX)
endif

BBLOCKDIR = $(CHPTDIR)/bblocks
BBLOCKNAME = bb$(NAMESUFFIX)

ifdef HASINFO
    FIELDINFO = M$(NM)V$(NV)
    POWERINFO = $(FIELDINFO)P$(NP)

    VERTEXDIR = $(CHPTDIR)/vertices
    VERTEXNAME = $(POWERINFO)$(NAMESUFFIX)
    VERTEXFILE = $(VERTEXDIR)/$(VERTEXNAME).hf

    ifdef BBLOCK
        BBLOCKFILE = $(BBLOCKDIR)/$(BBLOCKNAME)/$(BBLOCK)$(RANK)$(FIELDINFO).hf
    endif
endif
