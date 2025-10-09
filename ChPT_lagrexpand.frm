#-
* Generates vertex factors efficiently by inserted precomputed expansions of Lagrangian building blocks
* into the Lagrangian, using precomputed tables of integer partitions to get exactly the right orders
* with no waste or length computations, just a bit of system overhead.

off statistics;

#ifndef `ORDER'
    #message "ORDER (power-counting order) undefined"
    #terminate 1
#endif
#ifndef `NM'
    #message "NM (number of meson fields) undefined"
    #terminate 1
#endif
#ifndef `NV'
    #define NV "0"
#endif

#include- ChPT_bblocks.hf

functions <bb1>,...,<bb`ORDER'>;

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

* Load the Lagrangian
local vertex`SHORTNAME' =
    #include- ChPT_lagrangian.hf
    ;

#ifdef `SELECTLEC'
    if(match(`SELECTLEC')==0) discard;
#endif

.sort

* Hide the traces so that the Lagrangian is a pure prodoct of building blocks
#call untrace
#call hidetrace

* Prepend arguments to each building block indicating the power of each field
* This is done using tables of all ways to distribute N powers among NB building blocks
#procedure prependfield(FIELD)
    #do NBBLOCK=`ORDER',1,-1
        #system `MAKECMD' NBBLOCK=`NBBLOCK' `BBLOCKDIR'/rhs/`N`FIELD''on`NBBLOCK'.hf
        id,ifmatch->`FIELD'done <bb1?(?x1)>*...*<bb`NBBLOCK'?(?x`NBBLOCK')> =
            #include- `BBLOCKDIR'/rhs/`N`FIELD''on`NBBLOCK'.hf
            ;
    #enddo
    label `FIELD'done;
#endprocedure

#call prependfield(P)
#call prependfield(S)
#call prependfield(A)
#call prependfield(V)
#call prependfield(M)

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

*     Simplify QCD+QED case, repurpose A as the photon (not axial-vector) field
    #ifdef `QED'
        id V(?lz) = Q * A(?lz);
    #endif

*     Simplify equal-mass case
    #if ((isdefined(HASSCALAR)==0)&&(isdefined(FULLMASS)==0)&&(isdefined(PKEMASS)==0))
        id chi = m0p2;
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

*             Tracelessness
        id tr(f1?{<f1>,...,<f`NM'>}) = 0;
    #endif

*     Make the traces cyclic
*     This has little power since it isn't automatically combined with index
*      renaming, but it's better than nothing!
    id tr(?a) = Tr(?a);

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
    local unsub = vertex`SHORTNAME';
#endif

#call substituteall(vertex`SHORTNAME')
#call process(vertex`SHORTNAME')

#ifdef `CHECKTRANSFORM'
    #do TRANSFORM={CHIRAL,P,C}
        .sort
        local `TRANSFORM'vertex`SHORTNAME' = unsub;
        #call substituteall(`TRANSFORM'vertex`SHORTNAME')
        #call process(`TRANSFORM'vertex`SHORTNAME')
    #enddo
    .sort
    drop;

    #do TRANSFORM={CHIRAL,P,C}
        local [vertex`SHORTNAME':`TRANSFORM'violation] = `TRANSFORM'vertex`SHORTNAME' - vertex`SHORTNAME';
    #enddo

#endif

* * Use Cayley-Hamilton relations to simplify Nf=2 case
* #if (`NF'==2)&&(`NM'>0)
*     .sort
*
*     #define FLAVS "`NM'"
*     #include ChPT_external.hf
*
*     #call cayham(vertex`SHORTNAME')
*
*     sum <f1>,...,<f`NM'>;
* #endif

bracket i_,F,`LECS',Tr;
print +s;
.sort
#ifdef `CHECKTRANSFORM'
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
#else
    #write <`VERTEXFILE'> "* Generated by `NAME_' on `DATE_'"
    #write <`VERTEXFILE'> "%E",vertex`SHORTNAME'
#endif

.end

