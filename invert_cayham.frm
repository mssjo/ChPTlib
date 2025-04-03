* Cayley-Hamilton relations from Bijnens, Colangelo & Ecker '99
* With integrated tests
#-
off statistics;
#define NEXT "0"
#define PLAINCONV
#include- ChPTdiagram.hf

autodeclare symbol l,L, c,C,K, Y;


********************************************************************************
* 3 -> N

local lagrN = <K1*Y1>+...+<K115*Y115>;

#procedure cayham3
* SU(3) relations, eq. (B.1)
    id Y4 = 2*Y1 - 1/2*Y2 + Y3;
    id Y6 = 2/3*Y1 - 2/3*Y3 + Y5 + Y28 - 1/2*Y29
            - 1/3*Y30 - Y31 + 2/3*Y33 - 2/3*Y37 - 4/3*Y49
            + 1/3*Y50 - 2/3*Y52 + 8/3*Y54 - 1/3*Y57 + 2/3*Y58
            - 4/3*Y60 - 2/3*Y64 + 2/3*Y65 - 4/3*Y66 - 8/3*Y67
            + 2/3*Y68 + 4/3*Y86 - 2/3*Y87 - 4/3*Y88 + 2/3*Y89
            - 2/3*Y90 + 2/3*Y92 + 2*Y94 - 2/3*Y95 - 1/3*Y96
            - 2/3*Y97 + 1/3*Y99 + 2/3*Y105 + 1/3*Y106;
    id Y10 = -4*Y7 + 2*Y8 + 3*Y9 - 2*Y11 + 2*Y12;
    id Y15 = Y7 - 1/2*Y9 + Y11 - Y12 + Y13;
    id Y16 = 2*Y7 + Y8 - 3/2*Y9 + Y11 - Y12 + Y14;
    id Y22 = -4*Y19 + 4*Y20 + Y21 - 2*Y23 + 2*Y24;
    id Y32 = Y28 - 1/2*Y29 - Y30 + Y31;
    id Y36 = -4*Y33 + 4*Y34 + Y35 - 2*Y37 + 2*Y38;
    id Y51 = -4*Y49 + 5*Y50 - 2*Y52 + 2*Y53;
    id Y55 = Y49 - 1/2*Y50 + Y52 - Y53 + Y54;
    id Y56 = 2*Y49 - 1/2*Y50 + Y52 - Y53 + Y57;
    id Y59 = - Y49 + 3/4*Y50 - 3/2*Y52 + 3/2*Y53 + Y58;
    id Y61 = - Y49 + 1/4*Y50 + 1/2*Y52 + 3/2*Y53 -2*Y54 + 1/2*Y57 + Y60;
    id Y62 = - Y53 + 2*Y54 - 1/2*Y57 + Y60;
    id Y63 = Y49 - 3/4*Y50 + 3/2*Y52 - 5/2*Y53 + 4*Y54 - Y57 + Y60;
    id Y69 = Y64 - 1/2*Y65 + Y66;
    id Y70 = Y64 - Y65 + 2*Y67 + Y68;
    id Y74 = 2*Y71 - 1/2*Y72 + Y73;
    id Y79 = 2*Y75 + 2*Y76 - Y77 + Y78 - Y80;
    id Y93 = 2*Y90 - 1/2*Y91 + Y92;
    id Y98 = 2*Y94 + 2*Y95 - Y96 + Y97 - Y99;

#endprocedure
#procedure contact3
* eq. (3.4); Y116 is the SU(3) contact term
    id Y42 = 8*(Y116 - (1/12*Y25 - 1/8*Y26 + 1/24*Y27 + 1/4*Y39 - 1/8*Y40 - 1/4*Y41));
#endprocedure

#call cayham3
#call contact3

bracket Y1,...,Y116;
.sort
#do I=1,116
    local coeff`I' = lagrN[Y`I'];
#enddo
.sort

#define DROP3 "4,6,10,15,16,22,32,36,42,51,55,56,59,61,62,63,69,70,74,79,93,98"

* Construct SU(3) relations and the SU(3) lagrangian
local lagr3 =
    #define I3 "1"
    #define IN "1"
    #do DROP={`DROP3',117}
        #do I=`IN',`DROP'
            #if `I' != `DROP'
                + C`I3' * Y`IN'
                #write <invert_cayham_3toN.hf> "id C`I3' = %E;* Y`IN'",coeff`IN';
                #redefine I3 "{`I3'+1}"
            #endif
            #redefine IN "{`IN'+1}"
        #enddo
    #enddo
    ;

********************************************************************************
#ifdef `TEST3'

    .sort
    pushhide;

    local lagrNr = lagrN;
    #call renormalize(NNLO,2)

    id Nf = 3;
    id 1/Nf = 1/3;

    #do I=0,10
        id L`I'r = L`I';
    #enddo
    #do I=1,115
        id K`I'r = K`I';
    #enddo

    .sort
    pushhide;

    #undefine NFGENERAL
    #define NF "3"
    local lagr3r = lagr3;
    #call renormalize(NNLO,2)

    #do I=1,10
        id L`I'r = L`I';
    #enddo
    #do I=1,94
        id C`I'r = C`I';
    #enddo

    #call invertcayham(NNLO,3,N)
    #call invertcayham(NLO,3,N)

    .sort
    symbols a,b;
    pophide;
    drop lagr3r,lagrNr;

    local [should_be_zero] = lagr3r - lagrNr;

    bracket eps,F,Y1,...,Y116;
    print +s;
    .end
#endif
#ifdef `TEST3bis'
    .sort
    drop;

    local lagr = <K1*Y1>+...+<K115*Y115>;
    #call renormalize(NNLO,2)
    #call cayham3
    #call contact3
    id Nf=3;
    id 1/Nf = 1/3;

    if(count(eps,1)!=-2) discard;
    multiply -([d-4]/kappa*F)^2;
    #call expandD(2)

    bracket Y1,...,Y116, L0,...,L10;
    print;
    .end

#endif
********************************************************************************
* 2 -> 3 and 2 -> N

.sort

* Undo eq. (3.4)
id Y116 = 1/12*Y25 - 1/8*Y26 + 1/24*Y27 + 1/4*Y39 - 1/8*Y40 - 1/4*Y41 + 1/8*Y42;

.sort

#procedure cayham2()
* SU(2) relations, eq. (B.3)

    id Y2 = 2*Y1;
    id Y8 = 2*Y7;
    id Y9 = 2*Y7;
    id Y11 = Y7;
    id Y12 = 0;
    id Y14 = 2*Y13;
    id Y18 = 2*Y17;
    id Y21 = 2*Y19;
    id Y24 = Y19 - Y20 + Y23;
    id Y27 = -2*Y25 + 3*Y26;
    id Y29 = Y28;
    id Y30 = 0;
    id Y35 = 2*Y33;
    id Y38 = Y33 - Y34 + Y37;
    id Y42 = -2*Y39 + Y40 + 2*Y41;
    id Y45 = Y43 - Y44;
    id Y50 = 2*Y49;
    id Y52 = Y49;
    id Y53 = 0;
    id Y57 = 2*Y54;
    id Y60 = Y49 + Y54 - Y58;
    id Y64 = 2*Y67;
    id Y65 = 2*Y67;
    id Y68 = -Y66 - Y67;
    id Y72 = 2*Y71;
    id Y77 = Y75 + Y76;
    id Y80 = Y76 + 1/2*Y78;
    id Y82 = 2*Y81;
    id Y83 = 2*Y85;
    id Y84 = 2*Y85;
    id Y88 = 1/2*Y86 + 1/2*Y89;
    id Y91 = 2*Y90;
    id Y96 = Y94 + Y95;
    id Y99 = Y95 + 1/2*Y97;
    id Y103 = 2*Y102;
    id Y106 = 1/2*Y105;
    id Y108 = Y107;

*     Redo the SU(3) ones
    id Y4 = Y3 + Y1;
    id Y6 = YX;* This is the difficult one, involves EOM etc. so h -> f+, f-, chi- but not chi+
    id Y10 = 4*Y7;
    id Y15 = Y13 + Y7;
    id Y16 = 2*(Y13 + Y7);
    id Y22 = 2*Y20;
    id Y32 = Y31 + 1/2*Y28;
    id Y36 = 2*Y34;
    id Y42 = Y40 - 2*Y39 + 2*Y41;
    id Y51 = 4*Y49;
    id Y55 = Y54 + Y49;
    id Y56 = 2*(Y54 + Y49);
    id Y59 = Y58 + 2*Y54 + Y49;
    id Y61 = 3*Y49 + 2*Y54 - Y58;
    id Y62 = Y58 + Y54;
    id Y63 = Y58 + 2*Y54 + Y49;
    id Y69 = Y66 + Y67;
    id Y70 = 2*Y67;
    id Y74 = Y73 + Y71;
    id Y79 = (Y75 + Y76 + Y78)/2;
    id Y93 = Y92 + Y90;
    id Y98 = (Y94 + Y95 + Y97)/2;

#endprocedure
#procedure contact2
* eq. (3.5); Y117 is the SU(2) contact term
    id Y46 = 2*Y47 - Y48 - 4*Y113 + 2*Y117;
#endprocedure
#call cayham2
#call contact2

bracket Y1,...,Y117;
.sort
#do I=1,117
    local coeffN`I' = lagrN[Y`I'];
    local coeff3`I' = lagr3[Y`I'];
#enddo
.sort

#define DROP2 "2,4,6,8,...,12,14,...,16,18,21,22,24,27,29,30,32,35,36,38,42,45,46,50,...,53,55,56,57,59,...,65,68,69,70,72,74,77,79,80,82,83,84,88,91,93,96,98,99,103,106,108,116"

local lagr2 =
    #define I2 "1"
    #define IN "1"
    #do DROP={`DROP2',118}
        #do I=`IN',`DROP'
            #if `I' != `DROP'
                + c`I2' * Y`IN'
                #write <invert_cayham_2toN.hf> "id c`I2' = %E;* Y`IN'",coeffN`IN';
                #write <invert_cayham_2to3.hf> "id c`I2' = %E;* Y`IN'",coeff3`IN';
                #redefine I2 "{`I2'+1}"
            #endif
            #redefine IN "{`IN'+1}"
        #enddo
    #enddo
    ;

.sort

********************************************************************************
#ifdef `TEST2'

    .sort
    pushhide;

    local lagrNr = lagrN;

*     bracket Y1,...,Y117;
*     print +s;
*     .end
    #call renormalize(NNLO,2)

    id Nf = 2;
    id 1/Nf = 1/2;

    #do I=0,10
        id L`I'r = L`I';
    #enddo
    #do I=1,115
        id K`I'r = K`I';
    #enddo

    .sort
    pushhide;

    #undefine NFGENERAL
    #define NF "2"
    local lagr2r = lagr2;
    #call renormalize(NNLO,2)

    #do I=1,7
        id l`I'r = l`I';
    #enddo
    #do I=1,57
        id c`I'r = c`I';
    #enddo

    #call invertcayham(NNLO,2,N)
    #call invertcayham(NLO,2,N)

    .sort
    symbols a,b;
    pophide;
    drop;

    local [should_be_zero] = a*lagr2r - b*lagrNr;
*     Undo (3.5)
    id Y117 = (2*Y47 - Y48 - 4*Y113 - Y46)/2;

    id a = 1;
    id b = 1;
*     id Y?{Y1,...,Y117} = 1;

    if(count(eps,1)>-2) discard;

    bracket kappa,eps,F,Y1,...,Y116;
    print +s;
    .end
#endif
#ifdef `TEST2bis'
    .sort
    drop;

    local lagr = <K1*Y1>+...+<K115*Y115>;
    #call renormalize(NNLO,2)
    #call cayham2
    #call contact2
    id Nf=2;
    id 1/Nf = 1/2;

    if(count(eps,1)!=-2) discard;
    multiply -([d-4]/kappa*F)^2;
    #call expandD(2)

    bracket Y1,...,Y117, L0,...,L10;
    print;
    .end
#endif

print lagr2;
.end
