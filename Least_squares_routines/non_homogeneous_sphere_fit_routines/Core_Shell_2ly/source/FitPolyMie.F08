program MieScattFit

    use inputdata
    use MieScattMod
    use RI_funcMod
    use polyMod

    implicit none
    !############################################################################################
    !#######    Fit routine for a two-layer sphere.                                     #########
    !#######    Here we use the grid-search method to chi-squared minimization.         #########
    !#######    For more details, please read our paper:                                #########
    !#######    MIT License                                                             #########        
    !#######    Copyright (c) 2024 Christian Balderas, Rolando Castillo &               #########
    !#######    Complex-Fluids-IFUNAM                                                   #########
    !############################################################################################      
    REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
    integer :: i,j,k,status,&
            NRc,NFITVAR,ti,tf!,NANG
    real :: dummy1,dummy2,kappa,sigm,&
            Rc,Rcmax,Rcmin,Rcstep,stdgc,stdgcmax,stdgcmin,stdgcstep,Xc,&
            Rs,Rsmax,Rsmin,Rsminpar,Rsstep,stdgs,stdgsmax,stdgsmin,stdgsstep,Xs,&
            QSCA,chi2,QBACK,QEXT,diff1,dang,&
            Rcfit,Rsfit,stdgcfit,stdgsfit,kappafit,&
            trate
    COMPLEX :: nmax,nmin,nc,nstep,ncfit,nfit,&
            nsmax,nsmin,ns,nsstep,nsfit

    NAMELIST / PARFITCORE / nmax,nmin,nstep,Rcmax,Rcmin,Rcstep,stdgcmax,stdgcmin,stdgcstep
    NAMELIST / PARFITSHELL / Rsmax,Rsminpar,Rsstep,nsmax,nsmin,nsstep,stdgsmax,stdgsmin,stdgsstep
    
    Aly = 2
    ly = Aly

    !######################  Reading initial values and experimental data  #######################
    CALL inputs(WAVEL,REFMED,k0,HRAD,CHI2min)
    CALL ReadExpData(NANG,data)

    !##################################  Read fit conditions  ####################################
    OPEN(UNIT=1,FILE='parameters.nml',DELIM='APOSTROPHE')
    READ(UNIT=1,nml=PARFITCORE)
    READ(UNIT=1,nml=PARFITSHELL)
    CLOSE(UNIT=1)


    !##########################  Allocate memory for scattering arrays  ##########################
    CALL AllocScattArray(NANG)
    CALL AllocPolyArray(NANG)
    CALL AllocLyrsArrays(Aly)
    do i=1,NANG
        AANG(i) = rads(data(1,i))
    enddo
    
    !################################  Counting freedom's degrees  ###############################
    i = 0
    if (Rcmin .NE. Rcmax) then
        i = i+1
    endif
    if (nmax .NE. nmin) then
        i = i+1
    endif
    if (stdgcmax .NE. stdgcmin) then
        i = i+1
    endif
    if (Rsmin .NE. Rsmax) then
        i = i+1
    endif
    if (nsmax .NE. nsmin) then
        i = i+1
    endif
    if (stdgsmax .NE. stdgsmin) then
        i = i+1
    endif
    NFITVAR = i         !freedom's degrees          

    !########################################  HELLO WORLD!  ######################################
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',ACTION='WRITE')
    write(*,10)
    write(*,11) WAVEL,REFMED,k0,HRAD,k0*HRAD
    write(*,300)
    write(1,10)
    write(1,11) WAVEL,REFMED,k0,HRAD,k0*HRAD
    write(1,300)
    CLOSE(UNIT=1)
    
    !################################  START GRID-SEARCHING METHOD  ################################
    
    CALL system_clock(ti,trate)!Start chronometer
    Rc = Rcmin
    do while (Rc .LE. Rcmax) 
        Xc = k0*Rc      
        nc = nmax                
        do while (real(nc) .GE. real(nmin))
            REFREL = nc/REFMED          !Relative refractive index
            if (Rsminpar .EQ. 0.) then
                Rsmin = Rc + 1.    
                Rs = Rsmin
            else
                if (Rsminpar .GT. Rc) then
                    Rs = Rsminpar
                else
                    Rsmin = Rc + 1.    
                    Rs = Rsmin
                endif
            endif
            do while (Rs .LE. Rsmax)
                ns = nsmax
                do while(real(ns) .GE. real(nsmin))
                    if (stdgcmin .GT. 1.0 .AND. stdgcmin .EQ. stdgcmax) then
                        goto 1
                    endif
                    !######################  Multileyered-sphere monodisperse case  ######################
                    xl(1) = Rc
                    xl(2) = Rs
                    xl(:) = k0*xl(:)
                    ml(1) = nc
                    ml(2) = ns
                    ml(:) = ml(:)/REFMED
                    CALL MieMLy(Aly,ml,xl,QEXT,QSCA,QBACK,FF,NS1,NS2)
                    !write(*,*) QSCA,QEXT
                    kappa = 0.
                    do j=1,NANG
                        kappa = kappa + data(2,j)/NS1(j)
                    enddo
                    kappa = kappa/NANG
                    NS1(:) = kappa*NS1(:)
                    CALL chi2subrut(nang,nfitvar,data,NS1,chi2)
                    CALL system_clock(tf) !take time
                    !write(*,100) Rc,real(nc),1.00,Rs,real(ns),1.00,kappa,CHI2,(tf-ti)/trate
                    write(*,100,advance='no') Rc,real(nc),1.00,Rs,real(ns),1.00,kappa,CHI2,(tf-ti)/trate
                    write(*,'(A)',advance='no') '\r'
                    if (chi2 .LE. chi2min) then
                        CALL writechi2(chi2,Rc,nc,1.0,Rs,ns,1.0,kappa,&
                        chi2min,Rcfit,ncfit,stdgcfit,Rsfit,nsfit,stdgsfit,kappafit,(tf-ti)/trate)
                    endif
                    !#############################  Polydispersity over core  #############################
                    stdgc = stdgcmin
                    do while (stdgc .LE. stdgcmax)
                        dummy1 = Rc*stdgc**(-sqrt(18.))
                        dummy2 = Rc*stdgc**(sqrt(18.))
                        if (dummy1 .LE. 0. .OR. dummy2 .GE. Rs) then !Overlap flag
                            !write(*,*) "OVERLAP FLAG!!!"
                            exit
                        endif
                        CALL TwoLys_PolyRc(Rc,Rs,stdgc,QSCApoly,Fpoly,NS1poly,NS2poly)
                        kappa = 0.
                        do j=1,NANG
                            kappa = kappa + data(2,j)/NS1poly(j)
                        enddo
                        kappa= kappa/NANG
                        NS1poly(:) = kappa*NS1poly(:)
                        CALL chi2subrut(nang,nfitvar,data,Fpoly,chi2)
                        CALL system_clock(tf)
                        if (chi2 .LE. chi2min) then
                            CALL writechi2(chi2,Rc,nc,stdgc,Rs,ns,1.0,kappa,&
                            chi2min,Rcfit,ncfit,stdgcfit,Rsfit,nsfit,stdgsfit,kappafit,(tf-ti)/trate)
                        else
                            write(*,100,advance='no')  Rc,real(nc),stdgc,Rs,real(ns),1.00,kappa,CHI2,(tf-ti)/trate
                            write(*,'(A)',advance='no') '\r'
                            diff1 = abs(chi2 - chi2min)/chi2min
                            if (diff1 .GT. 0.5) then 
                                exit
                            endif
                        endif
                        stdgc = stdgc + stdgcstep
                    enddo
                    !#############################  Polydispersity over shell  #############################
                    xl(1) = Rc
                    xl(2) = Rs
                    xl(:) = k0*xl(:)
                    stdgs = stdgsmin
                    do while (stdgs .LE. stdgsmax)
                        dummy1 = Rs*stdgs**(-sqrt(18.))
                        if (dummy1 .LE. Rc) then !Overlap flag
                            write(*,*) 
                            write(*,*) "Overlap flag!"
                            exit
                        endif
                        CALL TwoLys_PolyRs(Rc,Rs,stdgs,QSCApoly,Fpoly,NS1poly,NS2poly)
                        kappa = 0.
                        do j=1,NANG
                            kappa = kappa + data(2,j)/NS1poly(j)
                        enddo
                        kappa= kappa/NANG
                        NS1poly(:) = kappa*NS1poly(:)
                        CALL chi2subrut(nang,nfitvar,data,NS1poly,chi2)
                        CALL system_clock(tf)
                        !write(*,100) Rc,real(nc),1.00,Rs,real(ns),stdgs,kappa,CHI2,(tf-ti)/trate
                        if (chi2 .LE. chi2min) then
                            CALL writechi2(chi2,Rc,nc,1.0,Rs,ns,stdgs,kappa,&
                            chi2min,Rcfit,ncfit,stdgcfit,Rsfit,nsfit,stdgsfit,kappafit,(tf-ti)/trate)
                        else
                            write(*,100,advance='no')  Rc,real(nc),1.0,Rs,real(ns),stdgs,kappa,CHI2,(tf-ti)/trate
                            write(*,'(A)',advance='no') '\r'
                            diff1 = abs(chi2 - chi2min)/chi2min
                            if (diff1 .GT. 0.5) then 
                                exit
                            endif
                        endif
                        stdgs = stdgs + stdgsstep
                    enddo
                    !###########################  Polydispersity over core and shell  ########################
1                   stdgc = stdgcmin
                    do while (stdgc .LE. stdgcmax)
                        stdgs = stdgsmin
                        do while(stdgs .LE. stdgsmax)    
                            dummy1 = Rs*stdgs**(-sqrt(18.))
                            dummy2 = Rc*stdgc**(sqrt(18.))
                            if (dummy1 .LE. dummy2) then !Overlap with shell's plateau
                                write(*,*)
                                write(*,*) "Overlap flag!"
                                exit
                            endif
                            dummy2 = Rc*stdgc**(-sqrt(18.)) !Overlap with core
                            if (dummy2 .LE. 0.) then
                                write(*,*)
                                write(*,*) "Overlap flag!"
                                exit
                            endif
                            xl(1) = Rc
                            xl(2) = Rs
                            xl(:) = k0*xl(:)
                            CALL TwoLys_PolyRcRs(Rc,Rs,stdgc,stdgs,nc,ns,QSCApoly,Fpoly,NS1poly,NS2poly)
                            !call sleep(1)
                            kappa = 0.
                            do j=1,NANG
                                kappa = kappa + data(2,j)/NS1poly(j)
                            enddo
                            kappa= kappa/NANG
                            NS1poly(:) = kappa*NS1poly(:)
                            CALL chi2subrut(nang,nfitvar,data,NS1poly,chi2)
                            CALL system_clock(tf)
                            if (chi2 .LE. chi2min) then
                                CALL writechi2(chi2,Rc,nc,stdgc,Rs,ns,stdgs,kappa,&
                                chi2min,Rcfit,ncfit,stdgcfit,Rsfit,nsfit,stdgsfit,kappafit,(tf-ti)/trate)
                            else
                                write(*,100,advance='no')  Rc,real(nc),stdgc,Rs,real(ns),stdgs,kappa,CHI2,(tf-ti)/trate
                                write(*,'(A)',advance='no') '\r'
                                diff1 = abs(chi2 - chi2min)/chi2min
                                if (diff1 .GT. 0.5) then 
                                    exit
                                endif
                            endif
                            stdgs = stdgs + stdgsstep
                        enddo
                        stdgc = stdgc + stdgcstep
                    enddo
                    ns = ns - nsstep
                enddo    
                Rs = Rs + Rsstep
            enddo 
            nc = nc - nstep
        enddo
        Rc = Rc + Rcstep
    enddo

    CALL system_clock(tf)

    !###################################  WRITING OUTPUT FILES  ###################################
    write(*,301)
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE') 
    write(1,301)
    CLOSE(UNIT=1)

    !recover fit values
    Rc = Rcfit
    nc = ncfit
    stdgc = stdgcfit
    Rs = Rsfit
    ns =  nsfit
    stdgs = stdgsfit
    kappa = kappafit

    NANG = 1800
    !Aly = 11
    CALL AllocScattArray(NANG)
    CALL AllocPolyArray(NANG)
    CALL AllocLyrsArrays(Aly)
    dang = 180./real(NANG)
    do i=1,NANG+1
        AANG(i) = rads(dang*real(i-1))
    enddo
    xl(1) = Rc
    xl(2) = Rs
    ml(1) = nc
    ml(2) = ns
    OPEN(UNIT=1,FILE='RIprofile_fitted.dat',STATUS='UNKNOWN',ACTION='WRITE')
    do i=1,Aly
        write(1,*) xl(i), real(ml(i))
    enddo
    CLOSE(UNIT=1)
    xl(:) = k0*xl(:)
    ml(:) = ml(:)/REFMED
    !CALL MieMLy(ly,ml,xl,QEXT,QSCA,QBACK,FF,NS1,NS2)
    if (stdgcfit .GT. 1.) then
        if (stdgsfit .EQ. 1.) then
            !WRITE(*,*) "POLYDISPERSITY ON Rc"
            CALL TwoLys_PolyRc(Rc,Rs,stdgc,QSCApoly,Fpoly,NS1poly,NS2poly)
        else
            !WRITE(*,*) "POLYDISPERSITY ON Rc AND Rs"
            CALL TwoLys_PolyRcRs(Rc,Rs,stdgc,stdgs,nc,ns,QSCApoly,Fpoly,NS1poly,NS2poly)
        endif
    else
        if (stdgsfit .GT. 1.) then
            CALL TwoLys_PolyRs(Rc,Rs,stdgs,QSCApoly,Fpoly,NS1poly,NS2poly)
            !WRITE(*,*) "POLYDISPERSITY ON Rs"
        else
            !WRITE(*,*) "MONODISPERSE PARAMETERS"
            CALL MieMLy(Aly,ml,xl,QEXT,QSCA,QBACK,FF,NS1,NS2)
            OPEN(UNIT=1,FILE='FF_fitted.dat',STATUS='UNKNOWN',ACTION='WRITE')
            OPEN(UNIT=2,FILE='S1_fitted.dat',STATUS='UNKNOWN',ACTION='WRITE')
            do i=1,NANG
                write(1,*) deg(AANG(i)),FF(i)
                write(2,*) deg(AANG(i)),NS1(i)
            enddo
            CLOSE(UNIT=2)
            CLOSE(UNIT=1)
            WRITE(*,302) Rc,real(nc),stdgc,Rs,real(ns),stdgs,CHI2min,QSCA,QSCA*pi*(rs**2.),kappa
            WRITE(*,*) "Total computing time = ",(tf-ti)/trate
            OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
            WRITE(1,302) Rc,real(nc),stdgc,Rs,real(ns),stdgs,CHI2min,QSCA,QSCA*pi*(rs**2.),kappa
            WRITE(1,*) "Total computing time = ",(tf-ti)/trate
            goto 2         
        endif
    endif
    OPEN(UNIT=1,FILE='FF_fitted.dat',STATUS='UNKNOWN',ACTION='WRITE')
    OPEN(UNIT=2,FILE='S1_fitted.dat',STATUS='UNKNOWN',ACTION='WRITE')
    do i=1,NANG
        write(1,*) deg(AANG(i)),Fpoly(i)
        write(2,*) deg(AANG(i)),NS1poly(i)
    enddo
    CLOSE(UNIT=2)
    CLOSE(UNIT=1)
    WRITE(*,302) Rc,real(nc),stdgc,Rs,real(ns),stdgs,CHI2min,QSCApoly,QSCApoly*pi*(rs**2.),kappa
    WRITE(*,*) "Total computing time = ",(tf-ti)/trate
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
    WRITE(1,302) Rc,real(nc),stdgc,Rs,real(ns),stdgs,CHI2min,QSCApoly,QSCApoly*pi*(rs**2.),kappa
    WRITE(1,*) "Total computing time = ",(tf-ti)/trate
2   CLOSE(UNIT=1)



10 FORMAT(  1x,"---------------------------------------------------------------------------------------------",/,&
            1x,"------------------------- Light scattering by single sphere program -------------------------",/,&
            1x,"---------------------------------------------------------------------------------------------",/) 
11 FORMAT(  1x,"---------------------------------------    INPUTS   -----------------------------------------",//,&
            10X,"WAVELENGTH (nm)=",f8.3,/,&
            10x,"REFMED = ",F8.4,/,&
            10X,"WAVE NUM, K=",F8.4,/,&
            10x,"HYDRODINAMICS RADIUS (nm)=",f8.3,/,&
            10X,"SIZE PARAMETER =",F8.4,//)
300 FORMAT( 1X,"-------------------------------- FITTING BY LEAST-SQUARED METHOD ------------------------------",/,&
            4X,'R_core',3X,'n_core',3x,'Poly_c',2X,'R_shell',3X,'n_shell',3x,'Poly_s',7x,'kappa',9x,'chi^2',&
                5X,'ctime(s)',/,&
            1x,"-----------------------------------------------------------------------------------------------")
100 FORMAT(2X,F6.1,5X,F5.3,5X,F4.2,2X,F6.1,5X,F5.3,6X,F4.2,3X,ES13.5,3x,F9.3,2X,F11.3)
101 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,F4.2,4X,F6.1,5X,F4.2,5X,'--Overlap--',3x,ES13.5)
102 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,F4.2,4X,F6.1,5X,F4.2,5X,'Error grows',3x,ES13.5)
103 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,F4.2,4X,F6.1,5X,F4.2,5X,'Too big Err',3x,ES13.5)
104 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,A,4X,F6.1,5X,A,3X,ES13.5,3x,ES13.5)
301 FORMAT(  1x,"--------------------------------- FITTING RESULTS ---------------------------------------------",/)
302 FORMAT( 10X,"Core Radius = ",F5.1,/,&
            10X,"Core Refractive Index = ",F5.3,/,&
            10x,"Core Polydispersity Stdgc =",F4.2,/,&
            10X,"Shell Radius = ",F5.1,/,&
            10x,"Shell Refractive Index = ",F5.3,/,&
            10x,"Shell Polydispersity Stdgs =",F4.2,/,&
            10x,"chi^2 fit parameter =",ES13.5,/,&
            10x,"Qsca =",ES13.5,/,&
            10x,"Csca (units^2) =", ES13.5,/,&
            10X,"I(theta) = Kappa*|S1|, kappa=", ES13.5,/,&
            /1x,"---------------------------------------- END --------------------------------------------------")

end program MieScattFit

subroutine chi2subrut(ndata,nfitvar,data,theory,chi2)
    implicit none
    !input variables
    integer,intent(in) :: ndata,nfitvar
    real,intent(in) :: data(3,ndata),theory(ndata)
    !output variables
    real,intent(out) :: chi2
    !local variables
    integer :: k

    real :: args

    chi2 = 0.
    do k=1,ndata
        chi2 = chi2 + ((data(2,k)-theory(k))/(data(3,k)))**2.
    enddo
    chi2 = chi2/real(ndata-Nfitvar)

    RETURN
end subroutine chi2subrut

subroutine writechi2(chi2,Rc,nc,stdgc,Rs,ns,stdgs,kappa,&
                    chi2min,Rcfit,ncfit,stdgcfit,Rsfit,nsfit,stdgsfit,kappafit,time)

    implicit none
    real,intent(in)     :: chi2,Rc,Rs,stdgc,stdgs,kappa,time
    complex,intent(in)  :: nc,ns
    real,intent(out)    :: chi2min,Rcfit,Rsfit,stdgcfit,stdgsfit,kappafit
    complex,intent(out) :: ncfit,nsfit

    chi2min = chi2
    Rcfit = Rc
    ncfit = nc
    stdgcfit =stdgc
    Rsfit = Rs
    nsfit = ns
    stdgsfit =stdgs
    kappafit = kappa
    !write(*,*)
    !write(*,*) nsfit, real(nsfit)
    !Rc,real(nc),1.00,Rs,real(ns),1.00,kappa,CHI2,(tf-ti)/trate
    WRITE(*,100) Rcfit,real(ncfit),stdgcfit,Rsfit,real(nsfit),stdgsfit,kappafit,CHI2,time
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
    WRITE(1,100) Rcfit,real(ncfit),stdgcfit,Rsfit,real(nsfit),stdgsfit,kappafit,CHI2,time
    CLOSE(UNIT=1)

100 FORMAT(2X,F6.1,5X,F5.3,5X,F4.2,2X,F6.1,5X,F5.3,6X,F4.2,3X,ES13.5,3x,F9.3,2X,F11.3)

    RETURN
end subroutine writechi2
