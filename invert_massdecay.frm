#-
off statistics;

#write <invert_massdecay.hf> "* Generated from `NAME_' on `DATE_'"
#write <invert_massdecay.hf> "* Follows the notation of Bijnens & Hermansson-Truedsson (2017)\n"

symbols M02,F0, M2,F, L0,L, m1,m2,m3,f1,f2,f3;
symbols n, kappa;
symbols x(:3), xi(:3);
autodeclare symbol aM,aF, bM,bF;

#define XLEXPAND(x,L,a) "(1
    #do I=1,3
        + `~x'^`~I' * (
        #do J=0,`~I'
            + `~L'^`~J' * `~a'`~I'`~J'
        #enddo
        )
    #enddo
    )"
#define NEXPAND(f,x) "(1
    + `~x'*(n*`~f'1)
    + `~x'^2*(n*`~f'2 + n*(n-1)/2 * `~f'1^2)
    + `~x'^3*(n*`~f'3 + n*(n-1)*`~f'1*`~f'2 + n*(n-1)*(n-2)/6 * `~f'1^3))"
#define LEXPAND(f,L,a) "
    #do I=1,3
        id `~f'`~I' =
        #do J=0,`~I'
            + `~L'^`~J' * `~a'`~I'`~J'
        #enddo
        ;
    #enddo
    "

local [M/M] = M02/M2 * `XLEXPAND(x,L0,aM)';
local [F/F] = F0/F * `XLEXPAND(x,L0,aF)';
local [x] = M02/F0^2 * kappa;


id M02^n? = M2^n * (`NEXPAND(m,xi)');
id F0^n?  = F^n * (`NEXPAND(f,xi)');
`LEXPAND(m,L,bM)'
`LEXPAND(f,L,bF)'

id M2/F^2 = xi/kappa;

bracket xi,L;
print +s [x];
.sort
drop [x];

id x = [x];

id L0 = L + xi*m1 + xi^2*(m2 - 1/2*m1^2) + xi^3*(m3 - m1*m2 + 1/3*m1^3);
`LEXPAND(m,L,bM)'

#do I=1,3
    bracket xi,L;
    .sort
    skip;

    #do J=0,`I'
        #do Q={M,F}
            local [b`Q'`I'`J'] = b`Q'`I'`J' - [`Q'/`Q'][xi^`I'*L^`J'];
        #enddo
    #enddo

    print;
    .sort
    skip; nskip [M/M],[F/F];

    #do J=0,`I'
        #do Q={M,F}
            #write <invert_massdecay.hf> "id b`Q'`I'`J' =%e",[b`Q'`I'`J'];
            id b`Q'`I'`J' = [b`Q'`I'`J'];
        #enddo
    #enddo
#enddo

.sort

local check = M2*([M/M]-1) + F*([F/F]-1);
.sort
#if termsin(check) != 0
    #message[massdecay]~~~ Inversion failed: the following should be 1
    skip; nskip [M/M],[F/F];
    print;
    .sort
    #remove invert_massdecay.hf
    #terminate
#endif
skip;
#include+ invert_massdecay.hf
.end
