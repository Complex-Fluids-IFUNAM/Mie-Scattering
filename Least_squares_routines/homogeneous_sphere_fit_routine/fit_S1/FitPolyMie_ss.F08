program MieScattFit

    use inputdata
    use MieScattMod
    use RI_funcMod
    use polyMod

    implicit none
    REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
    integer :: i,j,k,l,m,n,status,flag,flagout,flagchi1,&
            NRAD,NRc,NRs,Nstdgc,Nstdgs,NFITVAR!,NANG
    real :: dummy,dummy1,dummy2,kappa,Dns,RAD,Rcut,sigm,&
            Rc,Rcmax,Rcmin,Rcstep,stdgc,stdgcmax,stdgcmin,stdgcstep,Xc,&
            Rs,Rsmax,Rsmin,Rsstep,stdgs,stdgsmax,stdgsmin,stdgsstep,Xs,&
            QSCA,chi2,CHI2min,QBACK,QEXT,chi2loc,diff1,diff2,diffp,dang,&
            Rcfit,Rsfit,ncfit,stdgcfit,stdgsfit,kappafit
    COMPLEX :: nmax,nmin,nc,nstep,nfit,ifri1,ifri2

    NAMELIST / PARFIT / nmax,nmin,nstep,Rcmax,Rcmin,Rcstep,stdgcmax,stdgcmin,stdgcstep,chi2min

    CALL inputs(WAVEL,REFMED,k0,HRAD)
    CALL ReadExpData(NANG,data)

    !Allocate memory for scattering arrays
    !Aly = 500
    CALL AllocScattArray(NANG)
    CALL AllocPolyArray(NANG)
    !CALL AllocLyrsArrays(Aly)
    do i=1,NANG
        AANG(i) = rads(data(1,i))
    enddo
    OPEN(UNIT=1,FILE='inputs.nml',DELIM='APOSTROPHE')
    READ(UNIT=1,nml=PARFIT)
    CLOSE(UNIT=1)
    !## ENTER THE PARAMETERS FOR GRID-SEARCHING
    NFITVAR = 1         !Number of fit parameters
    !## refractive index rank
    !nmax = (1.59,0.)        
    !nmin = (1.59,0.)       
    !nstep = (0.01,0.)    
    !## Radii rank for the core
    !Rcmax = 150             
    !Rcmin = 10.            
    !Rcstep = 1.              
    NRc = (Rcmax-Rcmin)/Rcstep  
    !## Polydispersitive rank for the core
    !stdgcmax = 1.10              
    !stdgcmin = 1.01             
    !stdgcstep = 0.01            
    !############################################
    NStdgc = ((stdgcmax-stdgcmin)/stdgcstep)+1   
    !#### START GRID-SEARCHING METHOD
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',ACTION='WRITE')
    write(*,10)
    write(*,11) WAVEL,REFMED,k0,HRAD,k0*HRAD
    write(*,300)
    write(1,10)
    write(1,11) WAVEL,REFMED,k0,HRAD,k0*HRAD
    write(1,300)
    CLOSE(UNIT=1)
    !CHI2min = 1.E10      !Starting goodness of fit parameter "chi^2" 
    do i=1,NRc+1
        Rc = Rcmin + real(i-1)*Rcstep    !Set radius core
        X  = k0*Rc                      
        nc = nmax
        stdgc = 1.
        do while (real(nc) .GE. real(nmin))
            REFREL = nc/REFMED          !Relative refractive index
            CALL BHMIE(X,REFREL,NS1,NS2,FF,QSCA) !MONODISPERSE CASE
            kappa = 0.
            do j=1,NANG
                kappa = data(2,j)/NS1(j) + kappa
            enddo
            kappa = kappa/NANG
            NS1(:) = kappa*NS1(:)
            CALL chi2subrut(nang,nfitvar,data,NS1,chi2)
            if (chi2 .LE. chi2min) then
                chi2min = chi2
                Rcfit = Rc
                ncfit = nc
                stdgcfit = 1.
                kappafit = kappa
                WRITE(*,100) Rcfit,ncfit,stdgcfit,kappa,CHI2
                OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
                WRITE(1,100) Rcfit,ncfit,stdgcfit,kappa,CHI2
                CLOSE(UNIT=1)
            else
                write(*,100,advance='no') Rc,real(nc),stdgc,kappa,CHI2
                write(*,'(A)',advance='no') '\r' 
            endif
            stdgc = stdgcmin
            do while (stdgc .LE. stdgcmax) !Enter to polydispersity loop
                CALL polyscatt(stdgc,QSCApoly,Fpoly,NS1poly,NS2poly)
                kappa = 0.
                do j=1,NANG
                    kappa = data(2,j)/NS1poly(j) + kappa
                enddo
                kappa = kappa/NANG
                !kappa = data(2,1)/Fpoly(1)!0.5*(data(2,1)/Fpoly(1) + data(2,NANG)/Fpoly(NANG))
                NS1poly(:) = kappa*NS1poly(:)
                CALL chi2subrut(nang,nfitvar,data,NS1poly,chi2)
                if (chi2 .LE. chi2min) then
                    chi2min = chi2
                    Rcfit = Rc
                    ncfit = nc
                    stdgcfit = stdgc
                    kappafit = kappa
                    WRITE(*,100) Rcfit,ncfit,stdgcfit,kappa,CHI2
                    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
                    WRITE(1,100) Rcfit,ncfit,stdgcfit,kappa,CHI2
                    CLOSE(UNIT=1)
                else
                    write(*,100,advance='no') Rc,real(nc),stdgc,kappa,CHI2
                    write(*,'(A)',advance='no') '\r' 
                endif
                stdgc = stdgc + stdgcstep
            enddo
            nc = nc - nstep
        enddo
    enddo

    write(*,301)
    OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE') 
    write(1,301)
    CLOSE(UNIT=1)

    !writes output file
    Rc = Rcfit
    X = Rc*k0
    nc = ncfit
    REFREL = nc/REFMED
    stdgc = stdgcfit
!
    NANG = 1800
    CALL AllocScattArray(NANG)
    CALL AllocPolyArray(NANG)
    dang = 180./real(NANG)
    do i=1,NANG+1
        AANG(i) = rads(dang*real(i-1))
    enddo

    !################# WRITING FORM FACTOR FUNCTION TABLES ########################!

    if (stdgc .GT. 1.) then
        CALL polyscatt(stdgc,QSCApoly,Fpoly,NS1poly,NS2poly)
        OPEN(UNIT=1,FILE='FF_fitted_poly.dat',STATUS='UNKNOWN',ACTION='WRITE')
        OPEN(UNIT=2,FILE='S1_S2_fitted_poly.dat',STATUS='UNKNOWN',ACTION='WRITE')
            do i=1,NANG
                write(1,*) deg(AANG(i)),Fpoly(i)
                write(2,*) deg(AANG(i)),NS1poly(i),NS2poly(i)
            enddo
        CLOSE(UNIT=2)
        CLOSE(UNIT=1)
        !################ WRITING log's data ENDING ###########################################!
        WRITE(*,302) Rcfit,real(ncfit),stdgcfit,CHI2min,QSCApoly,QSCApoly*pi*(rc**2.),kappafit
        OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
        WRITE(1,302) Rcfit,real(ncfit),stdgcfit,CHI2min,QSCApoly,QSCApoly*pi*(rc**2.),kappafit
        CLOSE(UNIT=1)
    else
        OPEN(UNIT=1,FILE='FF_fitted_mono.dat',STATUS='UNKNOWN',ACTION='WRITE')
        OPEN(UNIT=2,FILE='S1_S2_fitted_mono.dat',STATUS='UNKNOWN',ACTION='WRITE')
        CALL BHMIE(X,REFREL,NS1,NS2,FF,QSCA)
        do i=1,NANG
            write(1,*) deg(AANG(i)),FF(i)
            write(2,*) deg(AANG(i)),NS1(i),NS2(i)
        enddo
        CLOSE(UNIT=2)
        CLOSE(UNIT=1)
        !################ WRITING log's data ENDING ###########################################!
        WRITE(*,302) Rcfit,real(ncfit),stdgcfit,CHI2min,QSCA,QSCA*pi*(rc**2.),kappafit
        OPEN(UNIT=1,FILE='log.out',STATUS='UNKNOWN',POSITION='APPEND',ACTION='WRITE')
        WRITE(1,302) Rcfit,real(ncfit),stdgcfit,CHI2min,QSCA,QSCA*pi*(rc**2.),kappafit
        CLOSE(UNIT=1)
    endif
    
    

10 FORMAT(  1x,"-----------------------------------------------------------",/,&
            1x,"-------- Light scattering by single sphere program --------",/,&
            1x,"-----------------------------------------------------------",/) 
11 FORMAT(  1x,"----------------------    INPUTS   ------------------------",//,&
            10X,"WAVELENGTH (nm)=",f8.3,/,&
            10x,"REFMED = ",F8.4,/,&
            10X,"WAVE NUM, K=",F8.4,/,&
            10x,"HYDRODINAMICS RADIUS (nm)=",f8.3,/,&
            10X,"SIZE PARAMETER =",F8.4,//)
300 FORMAT( 1X,"------------- FITTING BY LEAST-SQUARED METHOD -------------",/,&
            4x,'Radius',2X,'n_core',3X,'stdgc',8X,'kappa',11x,'chi^2',/,&
            1X,"-----------------------------------------------------------")
100 FORMAT(5X,F5.1,3x,F5.3,3X,F5.3,3x,ES13.5,3x,ES13.5)
101 FORMAT(5X,F5.1,3x,F5.3,3X,A)
301 FORMAT(  1x,"-----------------------------------------------------------",/)
302 FORMAT( 10X,"Fitted Radius(nm) = ",F5.1,/,&
            10x,"Fitted Refractive Index = ",F5.3,/,&
            10x,"Polydispersity Stdg =",F4.2,/,&
            10x,"chi^2 fit parameter =",ES13.5,/,&
            10x,"Qsca =",ES13.5,/,&
            10x,"Csca (nm^2) =", ES13.5,/,&
            10X,"I(theta) = |S1|*kappa, kappa=", ES13.5,/,&
            /1x,"-------------------------- END ----------------------------")

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
        chi2 = chi2 + ((data(2,k)-theory(k))/data(3,k))**2.
    enddo
    !chi2 = chi2/real(ndata-Nfitvar)

    RETURN
end subroutine chi2subrut