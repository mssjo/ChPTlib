#define MAKE "make --no-print-directory"

#procedure upart(TOT)
    .sort
    cfunction ipart;
    #if `TOT' <= 0
        global upart`TOT' = ipart();
    #else
        #do `HEAD'=1,{`TOT'-1}
            #if exists(upart{`TOT'-`HEAD'})==0
                #call upart({`TOT'-`HEAD'})
            #endif
            .sort
            skip;
            local upart`TOT'x`HEAD' =
                upart{`TOT'-`HEAD'}
                ;
            id ipart(?x) = ipart(`HEAD', ?x);
        #enddo
        .sort
        drop <upart`TOT'x1>,...,<upart`TOT'x{`TOT'-1}>;
        global upart`TOT' = <upart`TOT'x1>+...+<upart`TOT'x{`TOT'-1}>;
    #endif
#endprocedure

#procedure opart(TOT,MAXPART)
    .sort
    cfunction ipart;
    #if `TOT' <= 0
        global opart`TOT' = ipart();
    #else
        #if `MAXPART' < `TOT'
            #define MAX "`MAXPART'"
        #else
            #define MAX "`TOT'"
        #endif
        #do `HEAD'=1,{`MAX'-1},1
            #if exists(opart{`TOT'-`HEAD'}x`HEAD')==0
                #call opart({`TOT'-`HEAD'},`HEAD')
            #endif
            .sort
            skip;
            local opart`TOT'x`HEAD'x`MAX' =
                opart{`TOT'-`HEAD'}x`MAX'
                ;
            id ipart(?x) = ipart(`HEAD', ?x);
        #enddo
        .sort
        drop <opart`TOT'x1x`MAX'>,...,<opart`TOT'x{`TOT'-1}x`MAX'>;
        global opart`TOT'x`MAXPART' = <opart`TOT'x1x`MAX'>+...+<opart`TOT'x{`TOT'-1}x`MAX'>;
    #endif
#endprocedure

#procedure fpart(TOT,LEN)
    .sort
    cfunction ipart;
    #if `LEN' <= 0
        global fpart`TOT'x`LEN' = ipart();
    #elseif `LEN' == 1
        global fpart`TOT'x`LEN' = ipart(`TOT');
    #else
        #do `HEAD'=0,{`TOT'-1}
            #if exists(fpart{`TOT'-`HEAD'}x{`LEN'-1})==0
                #call fpart({`TOT'-`HEAD'},{`LEN'-1})
            #endif
            .sort
            skip;
            local fpart`TOT'x`HEAD'x{`LEN'-1} =
                fpart{`TOT'-`HEAD'}x{`LEN'-1}
                ;
            id ipart(?x) = ipart(`HEAD', ?x);
        #enddo
        .sort
        drop <fpart`TOT'x0x{`LEN'-1}>,...,<fpart`TOT'x{`TOT'-1}x{`LEN'-1}>;
        global fpart`TOT'x`MAXPART' = <fpart`TOT'x0x{`LEN'-1}>+...+<fpart`TOT'x{`TOT'-1}x{`LEN'-1}>;
    #endif
#endprocedure

#ifndef `TOT'
    #message[ipart]~~~ ERROR: Integer missing for ipart, should be called as 'form -d TOT=(integer) -d [OUF] ipart.frm
    #terminate
#elseif `TOT' < 0
    #message[ipart]~~~ ERROR: Integer must be non-negative (`TOT' provided)
    #terminate
#endif

#ifdef `O'
