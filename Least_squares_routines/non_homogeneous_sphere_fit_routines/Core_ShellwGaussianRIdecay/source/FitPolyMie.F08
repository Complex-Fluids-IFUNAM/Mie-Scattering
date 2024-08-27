program MieScattFit

    use inputdata
    use MieScattMod
    use RI_funcMod
    use polyMod

    implicit none
    !############################################################################################
    !#######    Fit routine for a multilayer sphere with a non-homogeneus refractive    #########
    !#######    index profile.                                                          #########
    !#######    The RI profile model is: Core-RI plateau-Gaussian decay                 #########
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
            QSCA,chi2,CHI2min,QBACK,QEXT,diff1,dang,&
            Rcfit,Rsfit,ncfit,stdgcfit,stdgsfit,kappafit,&
            trate
    COMPLEX :: nmax,nmin,nc,nstep,nfit

    NAMELIST / PARFITCOREOFSHELL / nmax,nmin,nstep,Rcmax,Rcmin,Rcstep,stdgcmax,stdgcmin,stdgcstep
    NAMELIST / PARFITSHELL / Rsmax,Rsminpar,Rsstep,stdgsmax,stdgsmin,stdgsstep
    NAMELIST / PARMLMIEFIT / Aly,CHI2min
 
    !Reading initial values and experimental data
    CALL inputs(WAVEL,REFMED,k0,HRAD,Rps,stdgps,nps)
    CALL ReadExpData(NANG,data)

    !Read parameters's intervals to fit
    OPEN(UNIT=1,FILE='parameters.nml',DELIM='APOSTROPHE')
    READ(UNIT=1,nml=PARFITCOREOFSHELL)
    READ(UNIT=1,nml=PARFITSHELL)
    READ(UNIT=1,nml=PARMLMIEFIT)
    CLOSE(UNIT=1)

    !Allocate memory for scattering arrays
    CALL AllocScattArray(NANG)
    CALL AllocPolyArray(NANG)
    CALL AllocLyrsArrays(Aly)
    do i=1,NANG
        AANG(i) = rads(data(1,i))
    enddo
    !## ENTER THE PARAMETERS FOR GRID-SEARCHING
    NFITVAR = 5.         !Number of fit parameters          

    !## HELLO WORLD!
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',ACTION='WRITE')
    write(*,10)
    write(*,11) WAVEL,REFMED,k0,HRAD,k0*HRAD
    write(*,300)
    write(1,10)
    write(1,11) WAVEL,REFMED,k0,HRAD,k0*HRAD
    write(1,300)
    CLOSE(UNIT=1)
    
    !#### START GRID-SEARCHING METHOD
    
    CALL system_clock(ti,trate)!Start chronometer
    Rc = Rcmin
    nc = nmax
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
    do while (Rc .LE. Rcmax) 
        Xc = k0*Rc                      
        do while (real(nc) .GE. real(nmin))
            REFREL = nc/REFMED          !Relative refractive index
            do while (Rs .LE. Rsmax)
                !#### Multileyered-sphere monodisperse case"
                sigm = sigmgauss(Rs,Rc,real(nc),REFMED+0.0001,REFMED)
                xl(1) = Rps
                ml(1) = nps
                CALL CSGaussProfv2(Rc,Rs,REFMED,nc,sigm,xl,ml,ly)
                xl(:) = k0*xl(:)
                ml(:) = ml(:)/REFMED
                CALL MieMLy(ly,ml,xl,QEXT,QSCA,QBACK,FF,NS1,NS2)
                kappa = 0.
                do j=1,NANG
                    kappa = kappa + data(2,j)/NS1(j)
                enddo
                kappa = kappa/NANG
                NS1(:) = kappa*NS1(:)
                CALL chi2subrut(nang,nfitvar,data,NS1,chi2)
                CALL system_clock(tf) !take time
                write(*,100,advance='no') Rps,1.0,Rc,real(nc),1.00,Rs,1.00,kappa,CHI2,(tf-ti)/trate
                write(*,'(A)',advance='no') '\r'
                if (chi2 .LE. chi2min) then
                    CALL writechi2(chi2,Rps,1.0,Rc,real(nc),Rs,1.0,1.0,kappa,&
                    chi2min,Rcfit,ncfit,Rsfit,stdgcfit,stdgsfit,kappafit,(tf-ti)/trate)
                endif
                !############## Polydispersity over shell's plateau ########!
                stdgc = stdgcmin
                do while (stdgc .LE. stdgcmax)
                    dummy1 = Rc*stdgc**(-sqrt(18.))
                    dummy2 = Rc*stdgc**(sqrt(18.))
                    if (dummy1 .LE. Rps .OR. dummy2 .GE. Rs) then !Overlap with core and Rs
                        exit
                    endif
                    CALL GaussSphPolyRc(Rc,Rs,nc,stdgc,QSCApoly,Fpoly,NS1poly,NS2poly)
                    kappa = 0.
                    do j=1,NANG
                        kappa = kappa + data(2,j)/NS1poly(j)
                    enddo
                    kappa= kappa/NANG
                    NS1poly(:) = kappa*NS1poly(:)
                    CALL chi2subrut(nang,nfitvar,data,Fpoly,chi2)
                    CALL system_clock(tf)
                    if (chi2 .LE. chi2min) then
                        CALL writechi2(chi2,Rps,1.0,Rc,real(nc),Rs,stdgc,1.0,kappa,&
                        chi2min,Rcfit,ncfit,Rsfit,stdgcfit,stdgsfit,kappafit,(tf-ti)/trate)
                    else
                        write(*,100,advance='no')  Rps,1.00,Rc,real(nc),stdgc,Rs,1.,kappa,CHI2,(tf-ti)/trate
                        write(*,'(A)',advance='no') '\r'
                        diff1 = abs(chi2 - chi2min)/chi2min
                        if (diff1 .GT. 1) then 
                            exit
                        endif
                    endif
                    stdgc = stdgc + stdgcstep
                enddo
                !############## Polydispersity over shell ########!
                stdgs = stdgsmin
                do while (stdgs .LE. stdgsmax)
                    dummy1 = Rs*stdgs**(-sqrt(18.))
                    if (dummy1 .LE. Rc) then !Overlap with PS-core
                        exit
                    endif
                    CALL GaussSphPolyRs(Rc,Rs,nc,stdgs,QSCApoly,Fpoly,NS1poly,NS2poly)
                    kappa = 0.
                    do j=1,NANG
                        kappa = kappa + data(2,j)/NS1poly(j)
                    enddo
                    kappa= kappa/NANG
                    NS1poly(:) = kappa*NS1poly(:)
                    CALL chi2subrut(nang,nfitvar,data,NS1poly,chi2)
                    CALL system_clock(tf)
                    if (chi2 .LE. chi2min) then
                        CALL writechi2(chi2,Rps,1.0,Rc,real(nc),Rs,1.0,stdgs,kappa,&
                        chi2min,Rcfit,ncfit,Rsfit,stdgcfit,stdgsfit,kappafit,(tf-ti)/trate)
                    else
                        diff1 = abs(chi2 - chi2min)/chi2min
                        if (diff1 .GT. 1.) then
                            exit
                        else
                           write(*,100,advance='no') Rps,1.00,Rc,real(nc),1.00,Rs,stdgs,kappa,CHI2,(tf-ti)/trate
                           write(*,'(A)',advance='no') '\r'
                        endif
                    endif
                    stdgs = stdgs + stdgsstep
                enddo
                !############## Polydispersity over shell's plateau and shell ########!
                stdgc = stdgcmin
                do while (stdgc .LE. stdgcmax)
                    stdgs = stdgsmin
                    do while(stdgs .LE. stdgsmax)    
                        dummy1 = Rs*stdgs**(-sqrt(18.))
                        dummy2 = Rc*stdgc**(sqrt(18.))
                        if (dummy1 .LE. dummy2) then !Overlap with shell's plateau
                            exit
                        endif
                        dummy2 = Rc*stdgc**(-sqrt(18.)) !Overlap with core
                        if (dummy2 .LE. Rps) then
                            exit
                        endif
                        CALL GaussSphPolyRcRs(Rc,Rs,nc,stdgc,stdgs,QSCApoly,Fpoly,NS1poly,NS2poly)
                        kappa = 0.
                        do j=1,NANG
                            kappa = kappa + data(2,j)/NS1poly(j)
                        enddo
                        kappa= kappa/NANG
                        NS1poly(:) = kappa*NS1poly(:)
                        CALL chi2subrut(nang,nfitvar,data,NS1poly,chi2)
                        CALL system_clock(tf)
                        if (chi2 .LE. chi2min) then
                            CALL writechi2(chi2,Rps,1.0,Rc,real(nc),Rs,stdgc,stdgs,kappa,&
                            chi2min,Rcfit,ncfit,Rsfit,stdgcfit,stdgsfit,kappafit,(tf-ti)/trate)
                        else
                            write(*,100,advance='no')  Rps,1.00,Rc,real(nc),stdgc,Rs,stdgs,kappa,CHI2,(tf-ti)/trate
                            write(*,'(A)',advance='no') '\r'
                            diff1 = abs(chi2 - chi2min)/chi2min
                            if (diff1 .GT. 1.) then
                                exit
                            endif
                        endif
                        stdgs = stdgs + stdgsstep
                    enddo
                    stdgc = stdgc + stdgcstep
                enddo
                Rs = Rs + Rsstep
            enddo
            Rs = Rsfit    
            nc = nc - nstep
        enddo
        nc = ncfit
        Rc = Rc + Rcstep
    enddo

    CALL system_clock(tf)

    CALL inputs(WAVEL,REFMED,k0,HRAD,Rps,stdgps,nps)
    write(*,301)
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE') 
    write(1,301)
    CLOSE(UNIT=1)

    !writes output file
    Rc = Rcfit
    !Xc = Rc*k0
    nc = ncfit
    Rs = Rsfit
    stdgc = stdgcfit
    stdgs = stdgsfit
    !chi2min = chi
    kappa = kappafit
!!
    NANG = 1800
    !Aly = 11
    CALL AllocScattArray(NANG)
    CALL AllocPolyArray(NANG)
    CALL AllocLyrsArrays(Aly)
    dang = 180./real(NANG)
    do i=1,NANG+1
        AANG(i) = rads(dang*real(i-1))
    enddo
    sigm = sigmgauss(Rs,Rc,real(nc),REFMED+0.0001,REFMED)
    CALL CSGaussProfv2(Rc,Rs,REFMED,nc,sigm,xl,ml,ly)
    xl(1) = Rps
    ml(1) = nps
    OPEN(UNIT=1,FILE='RIprofile_fitted.dat',STATUS='UNKNOWN',ACTION='WRITE')
    do i=1,ly
        write(1,*) xl(i), real(ml(i))
    enddo
    CLOSE(UNIT=1)
    xl(:) = k0*xl(:)
    ml(:) = ml(:)/REFMED
    if (stdgcfit .GT. 1.) then
        if (stdgsfit .EQ. 1.) then
            !WRITE(*,*) "POLYDISPERSITY ON Rc"
            CALL GaussSphPolyRc(Rc,Rs,nc,stdgc,QSCApoly,Fpoly,NS1poly,NS2poly)
        else
            !WRITE(*,*) "POLYDISPERSITY ON Rc AND Rs"
            CALL GaussSphPolyRcRs(Rc,Rs,nc,stdgc,stdgs,QSCApoly,Fpoly,NS1poly,NS2poly)
        endif
    else
        if (stdgsfit .GT. 1.) then
            CALL GaussSphPolyRs(Rc,Rs,nc,stdgs,QSCApoly,Fpoly,NS1poly,NS2poly)
            !WRITE(*,*) "POLYDISPERSITY ON Rs"
        else
            !WRITE(*,*) "MONODISPERSE PARAMETERS"
            CALL MieMLy(ly,ml,xl,QEXT,QSCA,QBACK,FF,NS1,NS2)
            OPEN(UNIT=1,FILE='FF_fitted.dat',STATUS='UNKNOWN',ACTION='WRITE')
            OPEN(UNIT=2,FILE='S1_fitted.dat',STATUS='UNKNOWN',ACTION='WRITE')
            do i=1,NANG
                write(1,*) deg(AANG(i)),FF(i)
                write(2,*) deg(AANG(i)),NS1(i)
            enddo
            CLOSE(UNIT=2)
            CLOSE(UNIT=1)
            WRITE(*,302) Rps,real(nps),stdgps,Rc,real(nc),stdgc,Rs,stdgs,CHI2min,QSCA,QSCA*pi*(rs**2.),kappa
            WRITE(*,*) "Total computing time = ",(tf-ti)/trate
            OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
            WRITE(1,302) Rps,real(nps),stdgps,Rc,real(nc),stdgc,Rs,stdgs,CHI2min,QSCA,QSCA*pi*(rs**2.),kappa
            WRITE(1,*) "Total computing time = ",(tf-ti)/trate
            goto 1
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
    WRITE(*,302) Rps,real(nps),stdgps,Rc,real(nc),stdgc,Rs,stdgs,CHI2min,QSCApoly,QSCApoly*pi*(rs**2.),kappa
    WRITE(*,*) "Total computing time = ",(tf-ti)/trate
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
    WRITE(1,302) Rps,real(nps),stdgps,Rc,real(nc),stdgc,Rs,stdgs,CHI2min,QSCApoly,QSCApoly*pi*(rs**2.),kappa
    WRITE(1,*) "Total computing time = ",(tf-ti)/trate
1   CLOSE(UNIT=1)



10 FORMAT(  1x,"-----------------------------------------------------------------------------------------------------",/,&
            1x,"------------------------- Light scattering by a multi-layer sphere program --------------------------",/,&
            1x,"-----------------------------------------------------------------------------------------------------",/) 
11 FORMAT(  1x,"-------------------------------------------    INPUTS   ---------------------------------------------",//,&
            10X,"WAVELENGTH (nm)=",f8.3,/,&
            10x,"REFMED = ",F8.4,/,&
            10X,"WAVE NUM, K=",F8.4,/,&
            10x,"HYDRODINAMICS RADIUS (nm)=",f8.3,/,&
            10X,"SIZE PARAMETER =",F8.4,//)
300 FORMAT( 1X,"---------------------------------- FITTING BY LEAST-SQUARED METHOD ----------------------------------",/,&
            4X,'R_ps',5X,'n_psc',3x,'R_cs',2X,'RI_core',3X,'Poly_c',3x,'R_shell',3x,'Poly_s',7x,'kappa',10x,'chi^2',&
                4X,'ctime(s)',/,&
            1x,"-----------------------------------------------------------------------------------------------------")
100 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,F4.2,4X,F6.1,5X,F4.2,3X,ES13.5,3x,F9.3,1X,F11.3)
101 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,F4.2,4X,F6.1,5X,F4.2,5X,'--Overlap--',3x,ES13.5)
102 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,F4.2,4X,F6.1,5X,F4.2,5X,'Error grows',3x,ES13.5)
103 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,F4.2,4X,F6.1,5X,F4.2,5X,'Too big Err',3x,ES13.5)
104 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,A,4X,F6.1,5X,A,3X,ES13.5,3x,ES13.5)
301 FORMAT(  1x,"----------------------------------------- FITTING RESULTS ------------------------------------------",/)
302 FORMAT( 10X,"Core Radius (nm) = ",F5.1,/,&
            10X,"Core Refractive Index = ",F5.3,/,&
            10x,"Core Polydispersity Stdgps =",F4.2,/,&
            10X,"Core of Shell Radius(nm) = ",F5.1,/,&
            10x,"Core of Shell Refractive Index = ",F5.3,/,&
            10x,"Core of Shell Polydispersity Stdgc =",F4.2,/,&
            10X,"Shell Radius(nm) = ",F5.1,/,&
            10x,"Shell Polydispersity Stdgc =",F4.2,/,&
            10x,"chi^2 fit parameter =",ES13.5,/,&
            10x,"Qsca =",ES13.5,/,&
            10x,"Csca (nm^2) =", ES13.5,/,&
            10X,"I(theta) = Kappa*|S1|, kappa=", ES13.5,/,&
            /1x,"----------------------------------------------- END -----------------------------------------------")

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

subroutine writechi2(chi2,Rps,stdgRps,Rc,nc,Rs,stdgc,stdgs,kappa,&
                    chi2min,Rcfit,ncfit,Rsfit,stdgcfit,stdgsfit,kappafit,time)
    implicit none
    real,intent(in)     :: chi2,Rps,stdgRps,Rc,nc,Rs,stdgc,stdgs,kappa,time
    real,intent(out)    :: chi2min,Rcfit,ncfit,Rsfit,stdgcfit,stdgsfit,kappafit
    
    chi2min = chi2
    Rcfit = Rc
    ncfit = nc
    Rsfit = Rs
    stdgcfit =stdgc
    stdgsfit =stdgs
    kappafit = kappa
    WRITE(*,100) Rps,stdgRps,Rc,real(nc),stdgc,Rs,stdgs,kappa,CHI2,time
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
    WRITE(1,100) Rps,stdgRps,Rc,real(nc),stdgc,Rs,stdgs,kappa,CHI2,time
    CLOSE(UNIT=1)

    100 FORMAT(2X,F6.1,5X,F4.2,2X,F6.1,3x,F5.3,5X,F4.2,4X,F6.1,5X,F4.2,3X,ES13.5,3x,F9.3,1X,F11.3)

    RETURN
end subroutine writechi2
