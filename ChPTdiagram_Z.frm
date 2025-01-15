#define NEXT "2"
#define LOOPMOMENTA "ell"
#include- ChPTdiagram.hf

local assignLO = assign(f1,f2);
#call makevertex(LO,2,LO)
local assignNLO = assign(f1,f2);
#call makevertex(NLO,2,NLO)
local assign4 = assign(f1,f2,iprop(ell),oprop(ell));
#call makevertex(4,4,LO)

bracket tr,F,i_;
print +s;
.sort
drop;

local M2contact = vertLO;
local M2counter = vertNLO;
local M2tadpole = int(ell) * vert4/2;

#call connectprops(M2contact,0)
#call connectprops(M2counter,0)
#call connectprops(M2tadpole,1)

#call reduceoneloop(M2tadpole,ell)

id p2 = -p1;

bracket tr,F,i_;
print +s;
.sort
drop;

local M2LO = M2contact;
local M2NLO = M2counter + M2tadpole;

bracket p1;
.sort
drop;

global ZLO = -i_ * M2LO[p1.p1];
global ZNLO = -i_ * M2NLO[p1.p1];

#call renormalize(ZLO);
#call renormalize(ZNLO);
#call massdecay(,LO)

multiply replace_(f1,i, f2,j);

bracket Tr,d_,Fpi,i_;
print +s;
.sort


#define NAME ""
#ifdef `NF'
    #redefine NAME "`NAME'_NF`NF'"
#endif
#ifdef `KALPAR'
    #redefine NAME "`NAME'_KALPAR"
#endif
#ifdef `GENPAR'
    #redefine NAME "`NAME'_GENPAR"
#endif

#write <extleg/ZLO`NAME'.hf> "%e",ZLO
#write <extleg/ZNLO`NAME'.hf> "%e",ZNLO
.end
