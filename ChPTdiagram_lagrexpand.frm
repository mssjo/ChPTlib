* Generates vertex factors efficiently by inserted precomputed expansions of Lagrangian building blocks
* into the Lagrangian, using precomputed tables of integer partitions to get exactly the right orders
* with no waste or length computations, just a bit of system overhead.

Off statistics;

#-
#ifndef `NP'
    #message "NP (power-counting order) undefined"
    #terminate 1
#endif
#ifndef `NM'
    #message "NM (number of meson fields) undefined"
    #terminate 1
#endif
#ifndef `NV'
    #define NV "0"
#endif

#include- ChPTdiagram_bblocks.hf

functions <bb1>,...,<bb`NP'>;

#ifdef `DRYRUN'
*     With DRYRUN defined, substitute just expands the building blocks but doesn't substitute them
    #do RANK=0,`BBLOCKMAXRANK'
        #do BBLOCK={`BBLOCKS`RANK''}
            #call substitute(`BBLOCK',`RANK')
        #enddo
    #enddo
    #message[bblock]~~~ DRYRUN: expanding fields but not substituting in Lagrangian
    .end
#endif

* Prepare tabulated partitions of expansion powers
#ifndef `NOMAKE'
    #do NB=1,`NP'
        #system `MAKECMD' NB=`NB' `BBLOCKDIR'/rhs/`NM'on`NB'.hf
        #system `MAKECMD' NB=`NB' `BBLOCKDIR'/rhs/`NV'on`NB'.hf
    #enddo
#endif

* Load the Lagrangian
local vertex`POWERINFO' =
    #include- ChPTdiagram_lagrangian.hf
    ;
#ifdef `NF'
    #if `NP'==4
        #if `NF'==3
            id L0 = 0;
        #elseif `NF'==2
            id L0 = 0;
            id L1 = l1/4;
            id L2 = l2/4;
            id L3 = 0;
            id L4 = l4/8;
            id L5 = 0;
            id L6 = (l3+l4)/16;
            id L7 = l7/-16;
            id L8 = 0;
            id L9 = l6/-2;
            id L10 = l5;
            id H1 = h1;
            id H2 = h2;
        #endif
    #endif
#endif

#ifdef `SELECTLEC'
    if(match(`SELECTLEC')==0) discard;
#endif

.sort

* Hide the traces so that the Lagrangian is a pure prodoct of building blocks
#call untrace
#call hidetrace

* Prepend two arguments to each building block, representing the number of meson and vector fields desired, respectively
* This is done using tables of all ways to distribute N powers among NB building blocks
#do NB=`NP',1,-1
    id,ifmatch->Mdone <bb1?(?x1)>*...*<bb`NB'?(?x`NB')> =
        #include- `BBLOCKDIR'/rhs/`NV'on`NB'.hf
        ;
#enddo
label Mdone;
#do NB=`NP',1,-1
    id,ifmatch->Vdone <bb1?(?x1)>*...*<bb`NB'?(?x`NB')> =
        #include- `BBLOCKDIR'/rhs/`NM'on`NB'.hf
        ;
#enddo
label Vdone;
#call unhidetrace


* Substitute all the building blocks
#procedure substituteall(EXPR)
    .sort
    skip; nskip `EXPR';
    #do RANK=0,`BBLOCKMAXRANK'
        #do BBLOCK={`BBLOCKS`RANK''}
            #call substitute(`BBLOCK',`RANK')
        #enddo
    #enddo
#endprocedure
#procedure process(EXPR)
    .sort
    skip; nskip `EXPR';

*     Apply some more stages of the transformations
    #ifdef `TRANSFORM'
        #if `TRANSFORM'==P
*         Apply parity transformation to all Lorentz indices: emu = epsilon(mu) = 1 if mu=0 else -1
            #do MU={`INDEXLIST'}
                repeat id bb?{`XBBLOCKS'}(?x, `MU', ?y) = bb(?x, e`MU', ?y) * e`MU';
                argument `XBBLOCKS';
                    id e`MU' = `MU';
                endargument;
            #enddo
        #elseif `TRANSFORM'==C
            id Tr(?x) = Tr(reverse_(?x));
            id i_ = -i_;* NOTE: all i's must be manifest at this point
        #endif
    #endif

*     Pull out photon field -- version-dependent!
    id v(?lz) = Q * A(?lz);

*     Simplify equal-mass case
    #if (isdefined(FULLMASS)==0)&&(isdefined(PKEMASS)==0)
        id chi = mp2;
    #endif

    id 1/sqrt2 = sqrt2/2;
    id sqrt2^2 = 2;



    #call trace

    #if `NM'>0
*         Pull out flavor indices
        cyclesymmetrize tr;
        #do I=`NM',1,-1
            id,once tr(?a, Phi(?lz), ?b) = tr(?a, f`I', ?b) * phi(f`I', ?lz);
        #enddo

*         Make the traces cyclic
*         This has little power since it isn't automatically combined with index renaming, but it's better than nothing!
        id tr(?a) = Tr(?a);

*             Tracelessness
        id Tr(f1?{<f1>,...,<f`NM'>}) = 0;
    #endif

*     Finalize the transformations
    #ifdef `TRANSFORM'
        .sort
        skip; nskip `EXPR';
        #if `TRANSFORM'==CHIRAL
            repeat;
                #do G={gL,gR,h}
                    id Tr(?a, `G'dag, `G', ?b) = Tr(?a, ?b);
                    id Tr(?a, `G', `G'dag, ?b) = Tr(?a, ?b);
                    repeat;
                        id Tr(?a, `G'dag, `G'(mu?), ?b) = -Tr(?a, `G'dag(mu), `G', ?b);
                        id Tr(?a, `G', `G'dag(mu?), ?b) = -Tr(?a, `G'(mu), `G'dag, ?b);
                    endrepeat;
                #enddo
            endrepeat;
        #elseif `TRANSFORM'==P
            #do MU={`INDEXLIST'}
                if((count(e`MU',1)==1)||(count(e`MU',1)>2));
                    print "ERROR: In %t";
                    exit "ERROR: Lorentz index `MU' present to erroneous power!";
                endif;
                id e`MU'^2 = 1;
            #enddo
        #endif
    #endif

    #ifdef `NF'
        id Tr() = `NF';
    #endif
#endprocedure

.sort
#ifdef `CHECKTRANSFORM'
    local unsub = vertex`POWERINFO';
#endif

#call substituteall(vertex`POWERINFO')
#call process(vertex`POWERINFO')

#ifdef `CHECKTRANSFORM'
    #do TRANSFORM={CHIRAL,P,C}
        .sort
        local `TRANSFORM'vertex`POWERINFO' = unsub;
        #call substituteall(`TRANSFORM'vertex`POWERINFO')
        #call process(`TRANSFORM'vertex`POWERINFO')
    #enddo
    .sort
    drop;

    #do TRANSFORM={CHIRAL,P,C}
        local [vertex`POWERINFO':`TRANSFORM'violation] = `TRANSFORM'vertex`POWERINFO' - vertex`POWERINFO';
    #enddo

#endif

bracket i_,F,`LECS',Tr;
print +s;
.sort

#ifndef `CHECKTRANSFORM'
    #write <`VERTEXFILE'> "%E",vertex`POWERINFO'
#else
    #do EXPR={`ACTIVEEXPRNAMES_'}
        #if termsin(`EXPR')
            #message[bblocks]~~~ tansformation check failed: `EXPR'
            #define FAIL
        #endif
    #enddo
    #ifdef `FAIL'
        #terminate
    #else
        #message[bblocks]~~~ tansformation check successful
    #endif
#endif

.end

