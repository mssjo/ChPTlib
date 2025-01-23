* Expands individual Lagrangian building blocks to the desired order in the fields.
* Intended as an auxiliary to ChPTdiagram_bblocks.frm

#-
#include- ChPTdiagram_bblocks.hf

* Load integer partitions
#procedure ipart(N,UO)
    #system `MAKECMD' partitions/`N'`UO'.hf
    #include- partitions/`N'`UO'.hf
#endprocedure
#procedure upart(N)
    #call ipart(`N',u)
#endprocedure
#procedure opart(N)
    #call ipart(`N',o)
#endprocedure

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
            #ifndef `GENPAR'
                #terminate[bbexpand]~~~ ERROR: nontrivial Phi-expansion without GENPAR
            #endif
            #if `NM'==0
                1;
            #elseif {`NM'%2}==0
                0;
            #else
                #call opart({`NM'-1})
                    ;
                #if `RANK' != 0
                    #message[bbexpand]~~~ ERROR: derivatives must be applied after expanding u -> Phi
                    #terminate
                #endif

*                 Prepend 1 + (number of 1's) to the partition, remove 1's
                id ipart(?x) = ipart(1, ?x);
                repeat id ipart(n?, ?x, 1) = ipart(n+1, ?x);
*                 Turn the partition into a product of traces,
*                  except the first which becomes a traceless product of untraced phi's
                id ipart(n?, ?x) = (Phi^n - 1/Nf*trPhi(n)) * ipart(?x) * a(n, ?x);
                repeat id ipart(?x, n?) = ipart(?x) * trPhi(n);
                id ipart() = 1;

                #ifndef `NFGENERAL'
                    id 1/Nf = 1/`NF';
                #endif

                id a(1) = 1;

                .sort:>>Phi GENPAR<<;

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
                #ifdef `KALPAR'
*                 Kalman's preferred parametrization
                    #if (`NM'==1)
                        i_ * `PM'Phi/(sqrt2*F);
                    #elseif ({`NM'%2}!=0)
                        0;
                    #else
                        (trPhi(2)/(4*F^2))^{`NM'/2} * kalmancoeff({`NM'/2});
                        id kalmancoeff(n?) = -binom_(2*n,n) / ( 2^(2*n) * (2*n-1) );

                        id trPhi(2) = tr(Phi,Phi);
                    #endif
                #else
*                 Exponential parametrization
                    (`PM'i_*Phi/(sqrt2*F))^`NM' * invfac_(`NM');

*                     ...with general reparametrization of the field
                    #ifdef `GENPAR'
                        id Phi^`NM' =
                            #call upart(`NM')
                            ;
                        repeat;
                            id ipart(n?odd_, ?a) = Phi(n,0)/F^(n-1) * ipart(?a);
                            al ipart(n?even_, ?a) = 0;
                        endrepeat;
                        id ipart() = 1;

                        #call substitute(Phi,0)

                    #endif
                #endif

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
            #ifdef `STRICT'
                #message[bbexpand]~~~ ERROR: unknown building block "`BBLOCK'"
                #terminate
            #else
                #message[bbexpand]~~~ WARNING: definition of "`BBLOCK'" not known, setting to zero
                0;
            #endif
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
