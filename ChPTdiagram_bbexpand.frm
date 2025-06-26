#-
* Expands Lagrangian building blocks to the desired order in the fields.
* Intended as an auxiliary to ChPTdiagram_lagrexpand.frm
#include- ChPTdiagram_bblocks.hf

#setexternalattr stderr=terminal

* Apply derivatives
#procedure deriv(MU)
    #define CONST "dR,trL,trR,chi,chidag"
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
* and transformation proerties given by X,Y:
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
        dL(mu) * `BBLOCK'(`NM',`NV'`INDICES') * dR(mu)
            #do Nm=0,`NM'
                #do Nv=0,`NV'
                    + Gamma(`Nm',`Nv', mu) * `BBLOCK'({`NM'-`Nm'},{`NV'-`Nv'}`INDICES')
                    - `BBLOCK'({`NM'-`Nm'},{`NV'-`Nv'}`INDICES') * Gamma(`Nm',`Nv', mu)
                #enddo
            #enddo
            ;
        #call substitute(`BBLOCK',`RANK')
        #call substitute(Gamma,1)
        #call deriv(mu)
    #else
        dL(mu) * `BBLOCK'(`NM',`NV'`INDICES') * dR(mu)
            #do Nm=0,`NV'
                - i_ * `X'(0,`Nm', mu) * `BBLOCK'(`NM',{`NV'-`Nm'}`INDICES')
                + i_ * `BBLOCK'(`NM',{`NV'-`Nm'}`INDICES') * `Y'(0,`Nm', mu)
            #enddo
            ;

        #if (`BBLOCK'==chi)||(`BBLOCK'==chidag)
            #call substitute`BBLOCK'
        #else
            #call substitute(`BBLOCK',`RANK')
        #endif
        #call substitute(lr,1)
        #call deriv(mu)
    #endif
#endprocedure

#procedure expanduGamma(PM)
    (
        #do Nx=0,`NM'
*             NOTE: we can't separate the NV=0 and NV=1 contributions since they mix under chiral transformations
            +    udag(`Nx',0,r) * (dL(mu) * u({`NM'-`Nx'},0,r) * dR(mu))
            `PM' u(`Nx',0,l) * (dL(mu) * udag({`NM'-`Nx'},0,l) * dR(mu))
            +    udag(`Nx',0,r) * ( -i_ * r(0,1, mu) * u({`NM'-`Nx'},0,r))
            `PM' u(`Nx',0,l) * ( -i_ * l(0,1, mu) * udag({`NM'-`Nx'},0,l))
        #enddo
        );

    #call substitute(lr,1)
    #call substitute(u,0)
    #call substitute(udag,0)
    #call deriv(mu)
#endprocedure

#procedure expandNGmatrix(PM)
* This handles the expansion of the Nambu-Goldstone matrix, which involves almost
* all dependence on PAR
    #if `PAR'==SQRT
*         The parametrization U = i_ Phi + sqrt(1 - Phi^2) (omitting normalization)
        #if (`NM'==1)
            i_ * `PM'Phi/(sqrt2*F);
        #elseif ({`NM'%2}!=0)
            0;
        #else
            (trPhi(2)/(4*F^2))^{`NM'/2} * coeff({`NM'/2});
            id coeff(n?) = -binom_(2*n,n) / ( 2^(2*n) * (2*n-1) );

            id trPhi(2) = tr(Phi,Phi);
        #endif
    #elseif `PAR'==GEN
*         The general parametrization is a bit special, since the field redefenition
*         causes mixing between different NM. We use partitions to help.

*         Rather than expanding exp(Phi) (omitting normalization) up to order NM
*          and then expanding the Phi's again, which would lead to lots of discarded terms,
*          list all ways of getting NM from up to NM Phi's, each expanded to just the
*          necessary order.
        #if `NM'>0
            #external ipart --form --odd --multiplicity `NM'
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
        (`PM'i_*Phi/(sqrt2*F))^`NM'
        #if `NM'>=2
            #switch `PAR'
                #case EXP
*                     Exponential parametrization
                    * invfac_(`NM');
                    #break
                #case CAY
*                     Cayley parametrization
                    / 2^{`NM'-1};
                    #break
                #case MIN
*                     Minimal parametrization
                    #if {`NM'%2}
                        * 0;
                    #else
                        * sign_({`NM'/2+1}) * / 2^{`NM'-2} * binom_({`NM'-2}, {`NM'/2-1})/`NM';
                    #endif
                    #break
                #case BIJ
*                     That parametrization Hans Bijnens used that one time
                    #if `NM'==2
                        / 2;
                    #elseif `NM'==3
                        / 8;
                    #else
                        #switch {`NM'%4}
                            #case 0
                            #case 2
                                * 0;
                                #break
                            #case 1
                                * -fac_({(`NM'-1)/2-2}) * invfac_({(`NM'-1)/4-1}) / 2^{7*((`NM'-1)/4) - 1};
                                #break
                            #case 3
                                * +fac_({(`NM'+1)/4-2}) / 2^{5*((`NM'+1)/4)-1};
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
* Most definitions are simply defined in terms of other building blocks, and these are loaded to
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
            #if `NM'==0
                1;
            #elseif {`NM'%2}==0
                0;
            #else
                #if `RANK' != 0
                    #message[bbexpand]~~~ ERROR: derivatives must be applied after expanding u -> Phi
                    #terminate
                #endif

*                 This implements footnote 11 in Bijnens, Husek & Sjo (2021)
*                 which is of course the best footnote ever
                #external ipart --form --ordered {`NM'-1}
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
            #if `RANK' == 0
                #do Nm=0,`NM'
                    #do Nv=0,`NV'
                        + u`dag'(`Nm',`Nv') * u`dag'({`NM'-`Nm'},{`NV'-`Nv'})
                    #enddo
                #enddo
                ;
                #call substitute(u`dag',0)
            #else
                #call covar(U`dag',{`RANK'-1},`rl',`lr')
            #endif
            #break

        #case udag
            #redefine PM "-"
*           INTENTIONAL FALL-THROUGH
        #case u
*             The basic chiral field u
            #if `RANK' == 0
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
                #if `NV'>0
                    0
                #else
                    #do Nm=0,`NM'
                        +    udag(`Nm',0,r) * chi(0,0)    * udag({`NM'-`Nm'},0,l)
                        `PM' u   (`Nm',0,l) * chidag(0,0) * u   ({`NM'-`Nm'},0,r)
                    #enddo
                #endif
                ;
                #call substitute(u,0)
                #call substitute(udag,0)
                #call substitute(chi,0)
                #call substitute(chidag,0)
            #endif
            #break

        #case fm
            #redefine PM "-"
*             INTENTIONAL FALL-THROUGH
        #case fp
            #if `RANK'>2
                #call covar(`BBLOCK',{`RANK'-1},h,h)
            #else
                #do Nx=0,`NM'
                    +    u(`Nx',0,l) * FL(0,`NV', mu,nu) * udag({`NM'-`Nx'},0,l)
                    `PM' udag(`Nx',0,r) * FR(0,`NV', mu,nu) * u({`NM'-`Nx'},0,r)
                #enddo
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
                #if `NM'!=0
                    0
                #else
                    dL(mu)*`lr'(0,1,nu)*dR(mu) - dL(nu)*`lr'(0,1,mu)*dR(nu)
                    - i_ * (`lr'(0,1,mu)*`lr'(0,1,nu) - `lr'(0,1,nu)*`lr'(0,1,mu))
                #endif
                ;
                #call substitute(lr,1)
                #call deriv(mu)
                #call deriv(nu)
            #endif
            #break

*         These only appear in the p6 lagrangian, no derivatives needed
        #case h
            u(`NM',`NV',mu,nu) + u(`NM',`NV',nu,mu);

            #call substitute(u,2)
            #break
        #case chimx
            #redefine pm "m"
            #redefine mp "p"
*             INTENTIONAL FALL-THROUGH
        #case chipx
            chi`pm'(`NM',`NV',mu) - i_/2 * (
                #do Nm=0,`NM'
                    #do Nm=0,`NV'
                        + chi`mp'(`Nm',`Nm')*u({`NM'-`Nm'},{`NV'-`Nm'},mu)
                        + u({`NM'-`Nm'},{`NV'-`Nm'},mu)*chi`mp'(`Nm',`Nm')
                    #enddo
                #enddo
                );

            #call substitute(chi`pm',1)
            #call substitute(chi`mp',0)
            #call substitute(u,1)
            #break

        #default
*             #ifdef `STRICT'
                #message[bbexpand]~~~ ERROR: unknown building block "`BBLOCK'"
                #terminate
*             #else
*                 #message[bbexpand]~~~ WARNING: definition of "`BBLOCK'" not known, setting to zero
*                 0;
*             #endif
    #endswitch
    ;

id 1/sqrt2 = sqrt2/2;
id sqrt2^2 = 2;

.sort
symbols Mtag,Vtag;

#call untrace

* Count number of fields
id Phi(?lz) = Phi(?lz) * Mtag;
id v(?lz) = v(?lz) * Vtag;

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
if((count(Mtag,1)!=`NM')||(count(Vtag,1)!=`NV'));
    print "[bbexpand] WARNING: inexact term detected (NM=`NM', NV=`NV'):";
    print "    %t";
    discard;
endif;

id Mtag = 1;
id Vtag = 1;

print +s;
.sort

* Write output file
#write <`BBLOCKFILE'> "%E",bb`BBLOCK'`RANK'`FIELDINFO'

.end
