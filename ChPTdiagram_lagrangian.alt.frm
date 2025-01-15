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

function tr;
cfunction Tr(cyclic);
function M,Q, t;
cfunction phi;

symbol F,B0, Nf, L0,...,L10;
symbol C1,...,C115;
symbol c1,...,c1862;

symbol sqrt2, invsqrt2;

function u, uu, du, Gamma;
symbol dag;
function dcov, dpart, dstop;
function chip, chim, chi;
function fp, fm, FL, FR;

* cbb = "chiral building blocks"
function bb, bbb, bbpart;
set cbb: u, chi, chip, chim
         fp, fm, FL, FR,
         bbpart
;
* fsm = "flavour space matrices"
set fsm: M,Q, chi, t;

index mu, nu, ro, si, mu1, mu2;
set lz: mu, nu, ro, si;

local vertex =
    #include ChPTdiagram_lagrangian.hf
    ;

* Place fields outside traces
* so that distributive addition can do its thing.
repeat id tr(?a, bb?cbb(?l)) = tr(?a) * bb(?l);

* Writes out fields in terms of their definitions.
* Up until this point, u(mu) is a separate entity from u.
* Here, we rewrite it so that u(mu) is just the mu-derivative of u.
* TODO in the future: general treatment of covariant derivatives!
id u(mu?) = i_  * (u(dag)*(u(mu) - i_*l(mu)*u) - u*(u(dag,mu) - i_*r(mu)*u(dag)));
id u(mu?, nu?) = i_ * (
                         u(dag,mu)*(u(nu) - i_*l(nu)*u) - u(mu)*(u(dag,nu) - i_*r(nu)*u(dag))
                       + u(dag)*(u(mu,nu) - i_*l(mu,nu)*u - i_*l(nu)*u(mu)) - u*(u(dag,mu,nu) - i_*r(mu,nu)*u(dag) - i_*r(nu)*u(mu))
                       + Gamma(mu) * (u(dag,mu)*(u(nu) - i_*l(nu)*u) - u(mu)*(u(dag,nu) - i_*r(nu)*u(dag)))
                       - (u(dag,mu)*(u(nu) - i_*l(nu)*u) - u(mu)*(u(dag,nu) - i_*r(nu)*u(dag))) * Gamma(mu)
                       );
id Gamma(mu?) = 1/2 * (u(dag)*(u(mu) - i_*l(mu)*u) + u*(u(dag,mu) - i_*r(mu)*u(dag)));

id chip = u(dag)*chi*u(dag) + u*chi(dag)*u;
id chim = u(dag)*chi*u(dag) - u*chi(dag)*u;

id fp(mu?,nu?) = u*FL(mu,nu)*u(dag) + u(dag)*FR(mu,nu)*u;
id fm(mu?,nu?) = u*FL(mu,nu)*u(dag) - u(dag)*FR(mu,nu)*u;

.sort
*#call print
*.end

* In the following, we will expand the Lagrangian into all terms that
* contain exactly NM instances of M, without first expanding every
* u, u(dag) to that order and discarding unwanted terms, a process which
* is very wasteful.


* Equips all chiral building blocks with a power counter
id bb?cbb(?a) = bb(?a, 0);

* Joins all traces for easier handling, but marks where the splits were.
repeat id tr(?a) * bb?cbb(?b) = tr(?a, bb(?b));
id tr(bb?cbb(?a), ?b) = tr(bb(split, ?a), ?b);
repeat id tr(?a) * tr(?b) = tr(?a, ?b);

* Equips the joined trace with a counter for the desired number of fields.
id tr(?a) = tr(`NM', ?a);

* Ensures that all u's with derivatives are expanded to at least 1st order.
repeat id tr(i?pos0_, ?a, u(?l, mu?lz, 0), ?b)
    = tr(i-1, ?a, u(?l, mu, -1), ?b);
* There wasn't enough powers to do that - term vanishes!
id tr(i?neg_, ?a) = 0;
* u(?l, i) is expanded to order i+(1 if derivative, 0 otherwise);
* shifting the zero simplifies things below.
* -1 was just a marker that the derivatives had "eaten" enough orders.
argument tr;
  id u(?l, -1) = u(?l, 0);
endargument;

.sort
* #call print
* .end

* Recursively enumerates all combinations of powers that sum to NM
* e.g. u(2)*u(0) + u(1)*u(1) + u(0)*u(2) for NM=2
id tr(i?, bb?cbb(?l, j?pos0_), ?a) = tr(bb(?l, j+i), ?a);
repeat;
  id tr(u(?l, i?pos_), bb?cbb(?ll, j?int_), ?a)
    = tr(u(?l, i-1), bb(?ll, j+1), ?a)
      + u(?l, i) * tr(bb(?ll, j), ?a);

  id tr(u(?l, 0), ?b) = u(?l, 0) * tr(?b);
  id tr(chi(?l, i?int_), bb?cbb(?ll, j?int_), ?a)
    = chi(?l, 0) * tr(bb(?ll, i+j), ?a);
endrepeat;
id tr(bb?cbb(?a)) = bb(?a);
id tr() = 1;

* Shifts the zeroes back
id u(?l, mu?lz, i?int_) = u(?l, mu, i+1);

* #call print
* .end

* Reinstates the traces
id bb?cbb(split, ?a) = tr() * bb(?a);

* Expands to exactly the right power.
* Notation: M(mu, nu, ..., i) = d_mu( d_nu ( ... ( M^i) ... ))
id u(dag, ?l, i?int_) = (-i_*invsqrt2/F)^i * M(?l, i) * param(i);
id u(     ?l, i?int_) = ( i_*invsqrt2/F)^i * M(?l, i) * param(i);
id chi(   ?l, i?pos_) = 0;
id chi(   ?l, 0     ) = chi(?l);

** This is the old, slow way - becomes awful at O(p^6)!
*id U(?l)    = sum_(i, 0, `NM', ( i_*sqrt2/F)^i*M(?l, i)*param(i));
*id Udag(?l) = sum_(i, 0, `NM', (-i_*sqrt2/F)^i*M(?l, i)*param(i));
*id M(mu?, ?l, 0) = 0;
*if(count(F,-1) != `NM')
*  discard;

id param(0) = 1;
id param(1) = 1;
#if(`PAR' == EXP)
** Exponential parametrisation
  id param(i?) = invfac_(i);
#elseif(`PAR' == CAY)
** Cayley parametrisation
  id param(i?) = 1/(2^(i-1));
#elseif(`PAR' == MIN)
** Minimal parametrisation
  id param(i?odd_) = 0;
  id param(i?even_) = sign_(1 + i/2)/(2^(i-2)) * binom_(i-2,i/2 - 1)/i;
#elseif(`PAR' == BIJ)
** Bijnens parametrisation
  id param(2) = 1/2;
  id param(3) = 1/8;

  id param(i?even_) = 0;
  id param(i?odd_) = param((i-1)/2);

  id param(i?even_) = -fac_(i-2)/(2^(7*i/2 - 1) * fac_(i/2-1));
  id param(i?odd_) = fac_((i+1)/2 - 2)/(2^(5*(i+1)/2 - 1));
#else
  id param(i?) = 1/0;
  #message Unsupported parametrisation!
  .end
#endif

*id sqrt2^2 = 2;
id invsqrt2^2 = 1/2;

* #call print
* .end
.sort
skip;
nskip `VERT';

id chi(dag) = chi;

* Puts stuff inside traces
repeat id tr(?a)*M?fsm(?b) = tr(?a,M(?b));

id tr(?a) = Tr(?a);

*#call print
*.end
.sort
off statistics;
skip;
nskip `VERT';

* Expands derivatives of M
* Repeatedly swaps between tr (easier to handle)
* and Tr (allows simplification by cyclicity)
#do i=1, {`NM'+`NP'}
  id Tr(?a) = tr(?a);
  id tr(?a, M(?l, 1), ?b) = tr(?a, M(?l), ?b);
  id tr(?a, M(mu?, j?{>1}), ?b) = tr(?a, M(mu), M(j-1), ?b)
                                + tr(?a, M, M(mu, j-1), ?b);
* Only second derivatives for now - need different treatment for generality
  id tr(?a, M(mu?, nu?, j?{>1}), ?b)
                                = tr(?a, M(mu, nu), M(j-1), ?b)
                                + tr(?a, M(mu), M(nu, j-1), ?b)
                                + tr(?a, M(nu), M(mu, j-1), ?b)
                                + tr(?a, M, M(mu, nu, j-1), ?b);
  id tr(?a, M(i?pos_), ?b) = tr(?a, M, M(i-1), ?b);
  id tr(?a, M(0), ?b) = tr(?a, ?b);
  id tr(?a) = Tr(?a);
  .sort
  skip;
  nskip `VERT';
#enddo

* #call print
* .end
.sort
on statistics;
skip;
nskip `VERT';

#message O(p^`NP') `NM'-point vertex factor:
#call print
.sort

#endprocedure

