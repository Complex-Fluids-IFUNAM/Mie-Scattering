program MieScatt

    use inputdata
    use MieScattMod

    implicit none
    !############################################################################################
    !#######  Mie Scattering program for monodisperse homogeneous or multilayer sphere  #########
    !#######  MIT License                                                               #########        
    !#######  Copyright (c) 2024 Christian Balderas, Rolando Castillo &                 #########
    !#######  Complex-Fluids-IFUNAM                                                     #########
    !############################################################################################      
    REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
    integer :: i,j,k,status,Aly
    real :: dang

    NAMELIST / HOMOSPHERE / Rp,np
    NAMELIST / MLSPHERE / Aly

    !#####  Reading initial values and experimental data  #####!
    CALL inputs(WAVEL,REFMED,k0,HOMOSPHE,NANG)

    !#####################  Hello world  ######################!
    OPEN(UNIT=2,FILE='log.out',STATUS='UNKNOWN',ACTION='WRITE')
    WRITE(2,10)
    WRITE(*,10)
    WRITE(2,11) WAVEL,REFMED,k0,HOMOSPHE,NANG
    WRITE(*,11) WAVEL,REFMED,k0,HOMOSPHE,NANG

    !#############  Allocating memory for arrays  #############!
    CALL AllocScattArray(NANG)
    dang = 180./real(NANG)
    do i=1,NANG+1
        AANG(i) = rads(dang*real(i-1))
    enddo

    !####################  MIE SCATTERING  ####################!
    if (HOMOSPHE .eqv. .true.) then !Call BHMie
        OPEN(UNIT=1,FILE='parameters.nml',DELIM='APOSTROPHE')
        READ(UNIT=1,nml=HOMOSPHERE)
        CLOSE(UNIT=1)   
        
        WRITE(*,12) Rp,Rp*K0,real(np),imag(np),real(np/REFMED),imag(np/REFMED)
        WRITE(2,12) Rp,Rp*K0,real(np),imag(np),real(np/REFMED),imag(np/REFMED)
        
        CALL BHMIE(Rp*k0,np/REFMED,NS1,NS2,FF,QSCA,QEXT,QBACK)

    else !Call Mie Scattering program for a multilayer sphere
        OPEN(UNIT=1,FILE='parameters.nml',DELIM='APOSTROPHE')
        READ(UNIT=1,nml=MLSPHERE)
        CLOSE(UNIT=1)

        CALL AllocLyrsArrays(Aly)
        CALL ReadRIprofile(Aly,xl,ml)
        WRITE(2,13) xl(Aly),xl(Aly)*k0,Aly
        WRITE(*,13) xl(Aly),xl(Aly)*k0,Aly
        
        xl(:) = xl(:)*k0
        ml(:) = ml(:)/REFMED
        CALL MieMLy(Aly,ml,xl,QEXT,QSCA,QBACK,FF,NS1,NS2)        

    endif

    !#######################  WRITING OUTPUT FILES  #######################!
    OPEN(UNIT=3,FILE='Ang_FormFactor.dat',STATUS='UNKNOWN',ACTION='WRITE')
    OPEN(UNIT=4,FILE='Ang_S1S1_S2S2.dat',STATUS='UNKNOWN',ACTION='WRITE')
    do i=1,NANG
        write(3,*) deg(AANG(i)),FF(i)
        write(4,*) deg(AANG(i)),NS1(i),NS2(i)
    enddo
    CLOSE(UNIT=3)
    CLOSE(UNIT=4)
    WRITE(2,14) QSCA,QEXT,QEXT-QSCA,QBACK
    WRITE(*,14) QSCA,QEXT,QEXT-QSCA,QBACK
    CLOSE(UNIT=2)

    !###############################  END  ################################!



!#################################################  FORMATS' LIST  #################################################!    
10 FORMAT(  1x,"---------------------------------------------------------------------------------------------",/,&
            1x,"------------------------- Light scattering by single sphere program -------------------------",/,&
            1x,"---------------------------------------------------------------------------------------------",/) 
11 FORMAT(  1x,"---------------------------------------    INPUTS   -----------------------------------------",//,&
            1X,"You set the units. Consider them in the results.",/,&
            1X,"WAVELENGTH =",f8.3,/,&
            1x,"REFMED = ",F8.4,/,&
            1X,"WAVE NUM, K=",F8.4,/,&
            1x,"IS A HOMOGENEOUS SPHERE? ",L2,/,&
            1X,"POINTS' NUMBER TO CALCULATE BETWEEN 0° AND 180° = ",I5,/)
12 FORMAT(  1x,"-------------------------------- HOMOGENEOUS SPHERE -----------------------------------------",/,&
            1X,"Particle's Radius = ",F5.1,/,&
            1x,"Size parameter =",F5.1,/,&
            1X,"Particle's Refractive Index = (",F5.3,",",F5.3,")",/,&
            1X,"Relative Refractive Index = (",F5.3,",",F5.3,")",/)
13 FORMAT(  1x,"-------------------------------- MULTILAYER SPHERE -----------------------------------------",/,&
            1X,"Particle's Radius = ",F5.1,/,&
            1x,"Size parameter =",F5.1,/,&
            1X,"Number of layers =",I5,/)
14 FORMAT(  1X,"----- MIE SCATTERING RESULTS -----",/,&
            1x,"Scattering efficency Qsca =",ES13.5,/,&
            1x,"Extintion efficency  Qext =",ES13.5,/,&
            1x,"Absortion efficency  Qabs =",ES13.5,/,&
            1x,"Backscattering eff. Qback =",ES13.5,/,&
            1X,"Form factor and scattering amplitudes as a function of scattering angle are in the .dat file",/,&
            1x,"--------------------------------------- END ------------------------------------------------")

end program MieScatt