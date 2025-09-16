#-
* Expands Lagrangian building blocks to the desired order in the fields.
* Intended as an auxiliary to ChPT_lagrexpand.frm
#include- ChPT_bblocks.hf

#setexternalattr stderr=terminal


#define CONST "dR,trL,trR"
#ifndef `HASSCALAR'
    #redefine CONST "`CONST',chi,chidag"
#endif

* Apply derivatives
#procedure deriv(MU)
    #call untrace
    repeat;
        id dL(`MU') * bb?!{`CONST'}(?lz) = bb(?lz) * dL(`MU') + bb(`MU',?lz);
        al dL(`MU') * dR(`MU') = 0;
        al dL(`MU') * chi?{`CONST'} = chi * dL(`MU');
    endrepeat;
    id dR(`MU') = 1;
    #call trace
#endprocedure

* Apply covariant derivative to building block BBLOCK with RANK Lorentz indices
* and transformation properties given by X,Y:
* either some combination of l,r for B -> gX B gYdag (D derivative) or h,h for B -> h B hdag (nabla derivative)
#procedure covar(BBLOCK,RANK,X,Y)

*   Generate the appropriate list of indices according to the rank, leaving mu for the derivative
    #define INDICES ""
    #define R "`RANK'"
    #do INDEX={`INDEXLISTNOMU'}
        #if `R'<=0
            #breakdo
        #endif
        #redefine INDICES "`INDICES',`INDEX'"
        #redefine R "{`R'-1}"
    #enddo

    #if (`X'==h)&&(`Y'==h)
        dL(mu) * `BBLOCK'(`FIELDS'`INDICES') * dR(mu)
            #define LOOPBODY "
                #if `~NSCALAR'
                #else
                    + Gamma(`~FIELDS', mu) * `BBLOCK'(`~FIELDSCOMPL'`INDICES')
                    - `BBLOCK'(`~FIELDSCOMPL'`INDICES') * Gamma(`~FIELDS', mu)
                #endif"
            #call fieldloop()
            ;

        #call substitute(`BBLOCK',`RANK')
        #call substitute(Gamma,1)
        #call deriv(mu)
    #else
        dL(mu) * `BBLOCK'(`FIELDS'`INDICES') * dR(mu)
            #define LOOPBODY "
                #if (`~NM')||(`~NSCALAR')
                #else
                    - i_ * `X'(`~FIELDS', mu) * `BBLOCK'(`~FIELDSCOMPL'`INDICES')
                    + i_ * `BBLOCK'(`~FIELDSCOMPL'`INDICES') * `Y'(`~FIELDS', mu)
                #endif"
            #call fieldloop()
            ;

        #call substitute(`BBLOCK',`RANK')
        #call substitute(lr,1)
        #call deriv(mu)
    #endif
#endprocedure

#procedure expanduGamma(PM)
    #show NNM
    #if `NNSCALAR'
        0
    #else
        #define LOOPBODY "
            #if `NNVECTOR'==0
                +    udag(`~MESONONLY',r) * (dL(mu) * u   (`~MESONCOMPL',r) * dR(mu))
                `PM' u   (`~MESONONLY',l) * (dL(mu) * udag(`~MESONCOMPL',l) * dR(mu))
            #elseif `NNVECTOR'==1
                +    udag(`~MESONONLY',r) * ( -i_ * r(`~VECTORONLY', mu) * u   (`~MESONCOMPL',r))
                `PM' u   (`~MESONONLY',l) * ( -i_ * l(`~VECTORONLY', mu) * udag(`~MESONCOMPL',l))
            #else
                + 0
            #endif"
        (
        #call mesonloop()
        )
    #endif
    ;

    #call substitute(lr,1)
    #call substitute(u,0)
    #call substitute(udag,0)
    #call deriv(mu)
#endprocedure

#procedure expandNGmatrix(PM)
* This handles the expansion of the Nambu-Goldstone matrix, which involves almost
* all dependence on PAR
    #if {`NNSCALAR'+`NNVECTOR'}
        0;
    #elseif `PAR'==SQRT
*         The parametrization U = i_ Phi + sqrt(1 - Phi^2) (omitting normalization)
        #if (`NNM'==1)
            i_ * `PM'Phi/(sqrt2*F);
        #elseif ({`NNM'%2}!=0)
            0;
        #else
            (trPhi(2)/(4*F^2))^{`NNM'/2} * coeff({`NNM'/2});
            id coeff(n?) = -binom_(2*n,n) / ( 2^(2*n) * (2*n-1) );

            id trPhi(2) = tr(Phi,Phi);
        #endif
    #elseif `PAR'==GEN
*         The general parametrization is a bit special, since the field redefenition
*         causes mixing between different NNM. We use partitions to help.

*         Rather than expanding exp(Phi) (omitting normalization) up to order NNM
*          and then expanding the Phi's again, which would lead to lots of discarded terms,
*          list all ways of getting NNM from up to NNM Phi's, each expanded to just the
*          necessary order.
        #if `NNM'>0
            #external ipart --form --odd --multiplicity `NNM'
            #fromexternal+
                ;
        #else
            ipart();
        #endif
        id ipart(?a) = ipart(?a) * (`PM'i_/(sqrt2*F))^nargs_(?a) * invfac_(nargs_(?a));
        repeat id ipart(n?, ?a) = Phi(n,0) * ipart(?a);
        id ipart() = 1;

        #call substitute(Phi,0)

    #else
*         All other parametrizations are given in terms of Taylor series
        (`PM'i_*Phi/(sqrt2*F))^`NNM'
        #if `NNM'>=2
            #switch `PAR'
                #case EXP
*                     Exponential parametrization
                    * invfac_(`NNM');
                    #break
                #case CAY
*                     Cayley parametrization
                    / 2^{`NNM'-1};
                    #break
                #case MIN
*                     Minimal parametrization
                    #if {`NNM'%2}
                        * 0;
                    #else
                        * sign_({`NNM'/2+1}) * / 2^{`NNM'-2} * binom_({`NNM'-2}, {`NNM'/2-1})/`NNM';
                    #endif
                    #break
                #case BIJ
*                     That parametrization Hans Bijnens used that one time
                    #if `NNM'==2
                        / 2;
                    #elseif `NNM'==3
                        / 8;
                    #else
                        #switch {`NNM'%4}
                            #case 0
                            #case 2
                                * 0;
                                #break
                            #case 1
                                * -fac_({(`NNM'-1)/2-2}) * invfac_({(`NNM'-1)/4-1}) / 2^{7*((`NNM'-1)/4) - 1};
                                #break
                            #case 3
                                * +fac_({(`NNM'+1)/4-2}) / 2^{5*((`NNM'+1)/4)-1};
                        #endswitch
                    #endif
                    #break
                #default
                    #message[bbexpand]~~~ ERROR: unrecognized parametrization PAR=`PAR'
                    #terminate
            #endswitch
        #else
            ;* coefficients 0 and 1 are always 1
        #endif
    #endif
#endprocedure


* Obtain the building block from its definition, taking some switchy shortcuts to avoid code duplication.
* Most building blocks are simply defined in terms of other building blocks, and these are loaded to
* exactly the right orders in the fields to smoothen the process.
#define PM "+"
#define pm "p"
#define mp "m"
#define lr "l"
#define rl "r"
#define dag ""
local bb`BBLOCK'`RANK'`FIELDINFO' =
    #switch `BBLOCK'
        #case Phi
*             Only needed in the general parametrization
            #if `PAR'!=GEN
                #message[bbexpand]~~~ ERROR: nontrivial Phi-expansion without PAR=GEN
                #terminate
            #endif
            #if `NNM'==0
                1;
            #elseif {`NNM'%2}==0
                0;
            #else
                #if `RANK' != 0
                    #message[bbexpand]~~~ ERROR: derivatives must be applied after expanding u -> Phi
                    #terminate
                #endif

*                 This implements footnote 11 in Bijnens, Husek & Sjo (2021)
*                 which is of course the best footnote ever
                #external ipart --form --ordered {`NNM'-1}
                #fromexternal+
                    ;
*                 Prepend 1 + (number of 1's) to the partition, remove 1's
                id ipart(?x) = ipart(1, ?x);
                repeat id ipart(n?, ?x, 1) = ipart(n+1, ?x);
*                 Turn the partition into a product of traces,
*                  except the first which becomes a traceless product of untraced phi's
                id ipart(n?, ?x) = (Phi^n - 1/`NF'*trPhi(n)) * ipart(?x) * a(n, ?x);
                repeat id ipart(?x, n?) = ipart(?x) * trPhi(n);
                id ipart() = 1;

                id a(1) = 1;

                .sort:>>Phi PAR=GEN<<;

*                 Turn trPhi into an explicit trace
                id trPhi(1) = 0;
                repeat id trPhi(n?pos_, ?x) = trPhi(n-1, Phi,?x);
                id trPhi(0, ?x) = tr(?x);

            #endif
            #break

        #case Udag
            #redefine dag "dag"
            #redefine lr "r"
            #redefine rl "l"
*           INTENTIONAL FALL-THROUGH
        #case U
            #if `NNSCALAR'
                0
            #elseif `RANK' == 0
                #if `NNVECTOR'
                    0
                #else
                    #define LOOPBODY "
                        + u`dag'(`~MESONONLY') * u`dag'(`~MESONCOMPL')"
                    #call fieldloop()
                    ;

                    #call substitute(u`dag',0)
                #endif
            #else
                #call covar(U`dag',{`RANK'-1},`rl',`lr')
            #endif
            #break

        #case udag
            #redefine PM "-"
*           INTENTIONAL FALL-THROUGH
        #case u
*             The basic chiral field u
            #if `NNSCALAR'
                0
            #elseif `RANK' == 0
                #call expandNGmatrix(`PM')
            #elseif `RANK'>1
                #call covar(u,{`RANK'-1},h,h)
            #else
*                 This is the main u(mu) field, simplified by using its similarity to Gamma
*                 No elegant fallthrough solution here, since FORM doesn't allow conditional #break
                i_ *
                #call expanduGamma(-)
            #endif
            #break
        #case Gamma
            1/2 *
            #call expanduGamma(+)
            #break

        #case chim
            #redefine PM "-"
*             INTENTIONAL FALL-THROUGH
        #case chip
            #if `RANK'>0
                #call covar(`BBLOCK',{`RANK'-1},h,h)
            #else
                #if `NNVECTOR'
                    0
                #else
                    #define LOOPBODY "
                        +    udag(`~MESONONLY',r) * chi   (`SCALARONLY') * udag(`~MESONCOMPL',l)
                        `PM' u   (`~MESONONLY',l) * chidag(`SCALARONLY') * u   (`~MESONCOMPL',r)"
                    #call mesonloop()
                    ;

                    #call substitute(u,0)
                    #call substitute(udag,0)
                    #call substitute(chi,0)
                    #call substitute(chidag,0)
                #endif
            #endif
            #break

        #case fm
            #redefine PM "-"
*             INTENTIONAL FALL-THROUGH
        #case fp
            #if `RANK'>2
                #call covar(`BBLOCK',{`RANK'-1},h,h)
            #else
                #define LOOPBODY "
                    +    u   (`~MESONONLY',l) * FL(`VECTORONLY', mu,nu) * udag(`~MESONCOMPL',l)
                    `PM' udag(`~MESONONLY',r) * FR(`VECTORONLY', mu,nu) * u   (`~MESONCOMPL',r)"
                #call mesonloop()
                ;
                #call substitute(u,0)
                #call substitute(udag,0)
                #call substitute(FL,2)
                #call substitute(FR,2)
            #endif
            #break

        #case FR
            #redefine lr "r"
*             INTENTIONAL FALL-THROUGH
        #case FL
            #if `RANK'>2
                #call covar(`BBLOCK',{`RANK'-1},`lr',`lr')
            #else
                #if (`NNM')||(`NNSCALAR')
                    0
                #elseif `NNVECTOR'==1
                    dL(mu)*`lr'(`VECTORONLY',nu)*dR(mu) - dL(nu)*`lr'(`VECTORONLY',mu)*dR(nu)
                #elseif `NNVECTOR'==2
                    - i_ * (
                        #define LOOPBODY "
                            + `lr'(`~VECTORONLY',mu)*`lr'(`~VECTORCOMPL',nu)
                            - `lr'(`~VECTORONLY',nu)*`lr'(`~VECTORCOMPL',mu)"
                        #call fieldloop()
                        )
                #else
                    0
                #endif
                ;
                #call substitute(lr,1)
                #call deriv(mu)
                #call deriv(nu)
            #endif
            #break

*         These only appear in the p6 lagrangian, no derivatives needed
        #case h
            #if `NNSCALAR'
                0;
            #else
                u(`FIELDS',mu,nu) + u(`FIELDS',nu,mu);
            #endif

            #call substitute(u,2)
            #break
        #case chimx
            #redefine pm "m"
            #redefine mp "p"
*             INTENTIONAL FALL-THROUGH
        #case chipx
            chi`pm'(`FIELDS',mu) - i_/2 * (
                #define LOOPBODY "
                    + chi`mp'(`~FIELDS')*u(`~FIELDSCOMPL',mu)
                    + u(`~FIELDS',mu)*chi`mp'(`~FIELDSCOMPL')"
                #call fieldloop
                );

            #call substitute(chi`pm',1)
            #call substitute(chi`mp',0)
            #call substitute(u,1)
            #break

*         chitil = (det(chi)/chi)^-1, used in some contact terms
        #case chitildag
            #redefine rl "l"
            #redefine lr "r"
            #redefine dag "dag"
*             INTENTIONAL FALL-THROUGH
        #case chitil
            #if `RANK'>0
                #call covar(`BBLOCK',{`RANK'-1},`rl',`lr')
            #elseif (`NNM')||(`NNVECTOR')
                0;
            #else
                #ifdef `CHIISSCALAR'
*                     chitilde simplifies to chi^(Nf-1) when it is scalar
                    (chi`dag'(`SCALARONLY'))^{`NF'-1};
                    #call substitute(chi`dag',0)
                #else
*                     otherwise, it is a primitive object
                    `BBLOCK'(`SCALARONLY');
                #endif
            #endif
            #break

        #default
            #message[bbexpand]~~~ ERROR: unknown building block "`BBLOCK'"
            #terminate
    #endswitch
    ;

id 1/sqrt2 = sqrt2/2;
id sqrt2^2 = 2;

.sort
symbols
    #do X={`FIELDTYPES'}
        #ifdef `HAS`X''
            ,`X'tag
        #endif
    #enddo
    ;

#call untrace

* Count number of fields
id Phi(?lz) = Phi(?lz) * Mtag;
#do X={`SPURIONS'}
    #ifdef `HAS`X''
        id `X'(?lz) = `X'(?lz) * `X'tag;
    #endif
#enddo

#ifdef `TRANSFORM'
    #if `TRANSFORM'==CHIRAL
*         Simplify transformations using unitarity
        repeat;
            #do G={gL,gR,h}
                id `G'dag*`G' = 1;
                id `G'*`G'dag = 1;
                repeat;
                    id `G'dag * `G'(mu?) = -`G'dag(mu) * `G';
                    id `G' * `G'dag(mu?) = -`G'(mu) * `G'dag;
                endrepeat;
            #enddo
        endrepeat;
    #endif
#endif


#call trace

.sort

* Filter number of fields
* NOTE: most building blocks get generated to exactly the right order, but u, Gamma, FL and FR
* get mixed powers of l,r and therefore require this step.
* It's also good as a sanity check.
if(0
    #do X={`FIELDTYPES'}
        #ifdef `HAS`X''
            ||(count(`X'tag,1)!=`N`X'')
        #endif
    #enddo
    );
    print "[bbexpand] ERROR: inexact term detected:";
    #do X={`FIELDTYPES'}
        #ifdef `HAs`X''
            print "[bbexpand]   (N`X'=`N`X'')";
        #endif
    #enddo
    print "    %t";
    exit "[bbexpand] Terminating due to inexact term.";
*     discard;
endif;

#do X={`FIELDTYPES'}
    #ifdef `HAS`X''
        id `X'tag = 1;
    #endif
#enddo

print +s;
.sort

* Write output file
#write <`BBLOCKFILE'> "%E",bb`BBLOCK'`RANK'`FIELDINFO'

.end
