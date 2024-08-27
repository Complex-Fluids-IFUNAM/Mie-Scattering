module MieScattMod
implicit none
!##############################################################################!
!#######                     MIE SCATTERING MODULE                      #######!
!####### This module contain the numerical method to calculate the      #######!
!####### radiative properties for a homogeneous and for amultilayered   #######!
!####### sphere.                                                        #######!
!####### MIT License                                                    #######!        
!####### Copyright (c) 2024 Christian Balderas, Rolando Castillo &      #######!
!####### Complex-Fluids-IFUNAM                                          #######!
!##############################################################################!      
!VARIABLES FOR MIE SCATTERING
REAL :: X
COMPLEX :: REFREL
REAL,ALLOCATABLE,DIMENSION(:) :: NS1,NS2,FF,AANG


CONTAINS
subroutine AllocScattArray(NANG)
    implicit none
    INTEGER,intent(IN) :: NANG

    !Allocate memory for size and RI arrays
    if (.NOT. ALLOCATED(AANG)) then
        ALLOCATE(AANG(NANG),NS1(NANG),NS2(NANG),FF(NANG))
    else
        DEALLOCATE(AANG,NS1,NS2,FF)
        ALLOCATE(AANG(NANG),NS1(NANG),NS2(NANG),FF(NANG))
    endif

    RETURN
end subroutine AllocScattArray

subroutine BHMIE(X,REFREL,NS1,NS2,FF,QSCA)
    use inputdata
    implicit none
    !################################################################################
    !###### BOHREN AND HUFFMAN'S ALGORITHM FOR MIE SCATTERING                       #
    !###### C.F. Bohren, D.R. Huffman,                                              #
    !###### Absorption and Scattering of Light by Small Particles,                  #
    !###### Wiley, Weinheim, Germany, 1998. https://doi.org/10.1002/9783527618156.  #
    !################################################################################
    !###### EN ESTA RUTINA SOLO CALCULO S1 Y S2 PARA OBTENER F(THETA)
    !###### PARA EL CASO DE LUZ DEPOLARIZADA
    !ARGUMENTOS
    !Input
    REAL, INTENT(IN) :: X
    COMPLEX, INTENT(IN) :: REFREL
    !Output
    REAL, INTENT(OUT) :: QSCA
    REAL, INTENT(OUT), DIMENSION(NANG) :: NS1,NS2,FF
    !Variables locales
    INTEGER :: i,j,jj,nmx
    REAL*8 :: CHI,CHI0,CHI1,DANG,DX,EN,FN,P,PII,PSI,PSI0,PSI1, &
              THETA,XSTOP,YMOD,DUMMY1,DUMMY2
    REAL*8, DIMENSION(NANG) :: AMU,PI,PI0,PI1,TAU
    COMPLEX(KIND=8) AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
    COMPLEX, DIMENSION(NANG) :: S1,S2
    COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: D

    !Constantes a double punto flotasnte
    PII = 4.D0*DATAN(1.D0)
    DX = X
    DREFRL=REFREL
    Y=X*DREFRL
    YMOD=ABS(Y)

    !Termino maximo de la expansion de coeficientes de Memoria
    XSTOP=X+4.*(X**0.3333)+2
    nmx = MAX(XSTOP,YMOD)+15

    !reservo memoria para D(:)
    ALLOCATE(D(nmx))

    !Inicializo funciones angulares
    !DANG = 0.5*PII/REAL(NANG0-1,8)
    do i=1,NANG
        THETA  = AANG(i) !TO DOUBLE FLOAT NUMBER
        AMU(i) = DCOS(THETA)
    enddo
    do i=1,NANG
        PI0(i)=0.
        PI1(i)=1.
    enddo
    !inicializo factores de forma S1 y S2
    do i=1,NANG
        S1(i)=(0.,0.)
        S2(i)=(0.,0.)
    enddo
    !Calculo derivada logaritmica con recurrencia descendente
    !iniciando con (0.,0.) en i=nmx
    D(nmx) = (0.,0.)
    do i=1,nmx-1
        EN = nmx-i+1
        D(nmx-i) = (EN/Y)-(1./(D(nmx-i+1)+EN/Y))
    enddo
    !Calculando las funciones de Bessel-Riccati con recurrencia
    !hacia arriba. Inicializo
    PSI0 = DCOS(DX)
    PSI1 = DSIN(DX)
    CHI0 = -SIN(X)
    CHI1 = COS(X)
    XI1 = DCMPLX(PSI1,-CHI1)
    P = -1.
    QSCA = 0.E0
    !GSCA = 0.E0
    !inicia el calculo de la serie para S1 y S2
    serie:do i=1,nmx
        EN = i
        FN = (2.E0*EN+1.)/(EN*(EN+1))
        PSI = (2.E0*EN-1.)*PSI1/DX - PSI0
        CHI = (2.E0*EN-1.)*CHI1/DX - CHI0
        XI = DCMPLX(PSI,-CHI)
        !Guarda valor anterior
        if (i > 1) then
            AN1 = AN
            BN1 = BN
        endif
        !Calcula coeficientes an y bn
        !En este caso PSI1(XI1) es el anterior a PSI(XI)
        AN = (D(i)/DREFRL + EN/DX)*PSI - PSI1
        AN = AN/((D(i)/DREFRL + EN/DX)*XI - XI1)
        BN = (DREFRL*D(i) + EN/DX)*PSI - PSI1
        BN = BN/((DREFRL*D(i) + EN/DX)*XI - XI1)
        QSCA = QSCA + (2.*EN + 1.)*(CDABS(AN)*CDABS(AN) + CDABS(BN)*CDABS(BN))
        !Calculo de amplitudes de forma DE 0 A 90 GRADOS
        do j=1,NANG
            PI(j) = PI1(j)
            TAU(j) = EN*AMU(j)*PI(j)-(EN+1.)*PI0(j)
            S1(j) = S1(j) + FN*(AN*PI(j) + BN*TAU(j))
            S2(j) = S2(j) + FN*(AN*TAU(j) + BN*PI(j))
        enddo
        !!Recorro elementos de la serie
        PSI0 = PSI1
        PSI1 = PSI
        CHI0 = CHI1
        CHI1 = CHI
        XI1 = DCMPLX(PSI1,-CHI1)
        !Calculo proximos elementos de pi(theta)
        do j=1,NANG
            PI1(j) = ((2.*EN+1.)*AMU(j)*PI(j) - (EN+1.)*PI0(j))/EN
            PI0(j) = PI(j)
        enddo
    enddo serie
    !OUTPUT
    QSCA = (2./(X*X))*QSCA
    do j=1,NANG
        NS1(j) = real(S1(j))**2. + imag(S1(j))**2.
        NS2(j) = real(S2(j))**2. + imag(S2(j))**2.
        FF(j)  = 0.5*(NS1(j) + NS2(j))
    enddo

    RETURN
end subroutine BHMIE

subroutine MieMLy(ly,ml,xl,QEXT,QSCA,QBACK,FF,NS1,NS2)
    use inputdata
    implicit none
    !####################################################################################
    !###### YANG'S ALGORITHM FOR MIE SCATTERING BY A MULTILAYER SPHERE                  #
    !###### Yang, W. (2003).                                                            #
    !###### Improved recursive algorithm for light scattering by a multilayered sphere  # 
    !###### Applied Optics, 42(9), 1710. https://doi.org/10.1364/AO.42.001710           #
    !###### Coded by Christian Balderas (2024)                                          #
    !###### MIT License                                                                 #        
    !###### Copyright (c) 2024 Christian Balderas, Rolando Castillo &                   #
    !###### Complex-Fluids-IFUNAM                                                       #
    !####################################################################################
    REAL*8, PARAMETER :: PII = 4.0d0*DATAN(1.d0) !Pi=3.1415
    !input variables
    INTEGER, INTENT(IN) :: ly
    COMPLEX, INTENT(IN) :: ml(ly)
    REAL, INTENT(IN) :: xl(ly)
    !output variables
    REAL,INTENT(OUT) :: QEXT,QSCA,QBACK
    REAL,INTENT(OUT), DIMENSION(NANG) :: FF,NS1,NS2
    !local variables
    INTEGER :: i,j,jj,l,k,NSTOP,NMAX,idummy1,idummy2
    REAL*8 :: P,THETA,DANG,RN,ANG,FN
    COMPLEX(kind=8) :: XBACK,AN,BN
    real*8, DIMENSION(ly) :: DXL
    COMPLEX(kind=8), DIMENSION(ly) :: DML
    REAL*8,DIMENSION(NANG) :: AMU,PI,PI0,PI1,TAU
    COMPLEX(kind=8), DIMENSION(NANG) :: S1,S2

    !## SERIE STOPPING NUMBER
    NSTOP = xl(ly) + 4.0*(xl(ly)**0.3333) + 2.
    !## NMAX Wiscombe criteria for downward recurrence of D1n(x)
    idummy1 = nint(sqrt(real(ml(1)*xl(1))**2. + imag(ml(1)*xl(1))**2.))
    NMAX = MAX(NSTOP,idummy1) 
    do l=2,ly
        idummy1 = nint(sqrt(real(ml(l)*xl(l))**2. + imag(ml(l)*xl(l))**2.))
        idummy2 = nint(sqrt(real(ml(l)*xl(l))**2. + imag(ml(l)*xl(l-1))**2.))
        NMAX = MAX(NMAX,idummy1,idummy2)
    enddo
    NMAX = NMAX + 15 
    !WRITE(*,100) NSTOP,NMAX

    !Passing arrays from floats to double floats 
    do i=1,ly
        DXL(i) = xl(i)
        DML(i) = ml(i)
    enddo
    !Starting array of cos(theta) ###!
    !DANG = 0.5*PII/REAL(NANG0-1,8)
    do i=1,NANG
        THETA  = AANG(i) !TO DOUBLE FLOAT NUMBER
        AMU(i) = DCOS(THETA)
    enddo
    do i=1,NANG
        PI0(i)=0.
        PI1(i)=1.
    enddo
    !inicializo factores de forma S1 y S2
    do i=1,NANG
        S1(i)=(0.,0.)
        S2(i)=(0.,0.)
    enddo
    !Starting series
    QSCA = 0.D0
    QEXT = 0.D0
    !QPR = 0.D0
    XBACK = DCMPLX(0.D0,0.D0)
    P = -1.
    !#### SERIE OVER AN AND BN COEFFICIENTS #####
    serie:do i=1,NSTOP
        RN = REAL(i,8)
        FN = (2.E0*RN+1.)/(RN*(RN+1))
        call AnBnCoeff(i,nmax,ly,dml,dxl,AN,BN)
        QSCA = QSCA + (2.*RN + 1.)*(CDABS(AN)*CDABS(AN)+CDABS(BN)*CDABS(BN))
        QEXT = QEXT + (2.*RN + 1.)*(REAL(AN,8)+REAL(BN,8))
        XBACK = XBACK + (2.*RN + 1.)*(-1.)**i *(AN-BN)
        !#### SERIE OVER AMPLITUD FUNCTIONS S1 AND S2 #####
        do j=1,NANG
            PI(j) = PI1(j)
            TAU(j) = RN*AMU(j)*PI(j)-(RN+1.)*PI0(j)
            S1(j) = S1(j) + FN*(AN*PI(j) + BN*TAU(j))
            S2(j) = S2(j) + FN*(AN*TAU(j) + BN*PI(j))
        enddo
        !RECURRENCE RELATION FOR PI(THETA)
        do j=1,NANG
            PI1(j) = ((2.*RN+1.)*AMU(j)*PI(j) - (RN+1.)*PI0(j))/RN
            PI0(j) = PI(j)
        enddo
    enddo serie
    !OUTPUT
    QSCA = (2./(xl(ly)*xl(ly)))*QSCA
    QEXT = (2./(xl(ly)*xl(ly)))*QEXT
    QBACK = XBACK*DCONJG(XBACK)
    QBACK = (1./(xl(ly)*xl(ly)))*QBACK
    do j=1,NANG
        NS1(j) = real(S1(j))**2. + imag(S1(j))**2.
        NS2(j) = real(S2(j))**2. + imag(S2(j))**2.
        FF(j) = (NS1(j)+NS2(j))*0.5d0    
    enddo

100 FORMAT(/,1X,"-------------------------------------------------",/,&
             1x,"--- WISCOMBE'S CRITERIA OF STOPING SERIES AND ---",/,&
             1X,"---         NMAX TO DOWNWARD RECURSION        ---",/,&
             1X,"---",10X,"Nstop =",I3,5X,"Nmax =",I3,9x,"---",/,&
             1X,"-------------------------------------------------")
    RETURN
end subroutine MieMLy

subroutine AnBnCoeff(n,nmax,ly,ml,xl,AN,BN)
    implicit none
    !input variables
    INTEGER,INTENT(IN) :: n,nmax,ly
    REAL*8, INTENT(IN), DIMENSION(ly) :: xl 
    COMPLEX(KIND=8), INTENT(IN), DIMENSION(ly) :: ml
    !output variables
    COMPLEX(kind=8), INTENT(OUT) :: AN,BN
    !locar variables
    integer :: i,j,k,l
    REAL*8 :: dn,xlm1,xly,xlym1,PI,PI0,PI1,TAU
    COMPLEX(kind=8) :: z,z1,z2,Hal,Hbl,Halm1,Hblm1,D1nl,D1nlm1,D3nl,D3nlm1,Qln,dmmy,&
                        Ga1,Ga2,Gb1,Gb2,PSIn,PSInm1,XIn,XInm1
    
    !#############################################
    !Recurrence is over l-layers
    !Starting upward recurrence 
    !with Ha_n(m_1,x_1)=Hb_n(m_1,x_1)=D1n(m_1,x_1)
    !#############################################
    z1=ml(1)*xl(1)      !l=1
    call D1nRecu(z1,n,nmax,D1nl,dmmy)
    Halm1 = D1nl
    Hblm1 = D1nl
    do l=2,ly
        z1=ml(l)*xl(l-1)
        z2=ml(l)*xl(l)
        call D1nRecu(z1,n,nmax,D1nlm1,dmmy) !D1_n(m_l*x_(l-1))
        call D3nRecu(z1,n,nmax,D3nlm1,dmmy) !D3_n(m_l*x_(l-1))
        Ga1 = ml(l)*Halm1 - ml(l-1)*D1nlm1
        Ga2 = ml(l)*Halm1 - ml(l-1)*D3nlm1
        Gb1 = ml(l-1)*Hblm1 - ml(l)*D1nlm1
        Gb2 = ml(l-1)*Hblm1 - ml(l)*D3nlm1
        xlym1=xl(l-1)
        xly=xl(l)
        call QRecu(z1,z2,xlym1,xly,n,nmax,Qln)
        call D1nRecu(z2,n,nmax,D1nl,dmmy) !D1_n(m_l*x_l)
        call D3nRecu(z2,n,nmax,D3nl,dmmy) !D3_n(m_l*x_l)
        Hal = (Ga2*D1nl - Qln*Ga1*D3nl)/(Ga2 - Qln*Ga1)
        Hbl = (Gb2*D1nl - Qln*Gb1*D3nl)/(Gb2 - Qln*Gb1)
        !WRITE(*,*) Hal,Hbl
        Halm1 = Hal
        Hblm1 = Hbl
    enddo
    !Evaluation of AN and BN coefficients
    z = COMPLEX(xl(ly),0.d0)
    call PsiRecu(z,n,nmax,PSIn,PSInm1)  !PSI_n(x_ly) !PSI_(n-1)(x_ly)
    call XiRecu(z,n,nmax,XIn,XInm1)   !XI_n(x_ly)  !XI_(n-1)(x_ly)
    dn = n
    an = (Hal/ml(ly) + dn/xl(ly))*PSIn - PSInm1
    an = an/((Hal/ml(ly) + dn/xl(ly))*XIn - XInm1)
    bn = (ml(ly)*Hbl + dn/xl(ly))*PSIn - PSInm1
    bn = bn/((ml(ly)*Hbl + dn/xl(ly))*XIn - XInm1)
    
    RETURN
end subroutine AnBnCoeff

subroutine D1nRecu(z,n,nmax,D1n,D1nm1)
    implicit none
    !input variables
    INTEGER, INTENT(IN) :: nmax,n
    COMPLEX(kind=8), INTENT(IN) :: z
    !output variables
    COMPLEX(kind=8), INTENT(OUT) :: D1n,D1nm1
    !Local variables
    INTEGER :: i,j,k,l,Dn
    COMPLEX(kind=8) :: D1i,D1im1

    !#### Downward recurrence for D1_n(z) Log-derivative of PSI(z)
    D1i = (0.d0,0.d0) !i = nmax+1
    Dn = nmax
    do i=1,nmax-n+1 
        D1im1 = D1i + (Dn/z)
        D1im1 = 1.d0/D1im1
        D1im1 = (Dn/z) - D1im1
        D1n = D1i
        D1nm1 = D1im1
        D1i = D1im1  !Recurrence      
        Dn = nmax-i
    enddo

    RETURN
end subroutine D1nRecu

subroutine D3nRecu(z,n,nmax,D3n,d3nm1)
    implicit none
    !input variables
    INTEGER, INTENT(IN) :: n,nmax
    COMPLEX(kind=8), INTENT(IN) :: z
    !output variables
    COMPLEX(kind=8), INTENT(OUT) :: D3n,D3nm1
    !Local variables
    INTEGER :: i,j,k,l,Dn
    REAL*8 :: a,b
    COMPLEX(kind=8) :: D3i,D3im1,PXim1,PXi,D1n,D1nm1
    
    !#### Upward recurrence for D3_n(z) Log-derivative of Xi(z)
    a = REAL(z)
    b = imag(z)
    D3im1 = (0.d0,1.d0) !i-1=0
    PXim1 = COMPLEX(cos(2*a),sin(2*a)) 
    PXim1 = 0.5*(1. - PXim1*exp(-2*b))
    do i=1,n
        call D1nRecu(z,i,nmax,D1n,D1nm1)
        !Psi_i(z)*Xi_i(z) recurrence function
        PXi = PXim1*(i/z - D1nm1)*(i/z - D3im1)
        !D3_n(z) recurrence function
        D3i = D1n + (0.d0,1.d0)/PXi
        D3n = D3i
        D3nm1 = D3im1
        D3im1 = D3i
        PXim1 = PXi
    enddo

    RETURN
end subroutine D3nRecu

subroutine PsiRecu(z,n,nmax,PSIn,PSInm1)
    implicit none
    !input variables
    INTEGER, INTENT(IN) :: n,nmax
    COMPLEX(kind=8), INTENT(IN) :: z
    !output variables
    COMPLEX(kind=8), INTENT(OUT) :: PSIn,PSInm1
    !Local variables
    INTEGER :: i,j,k,l,Dn
    REAL*8 :: a,b
    COMPLEX(kind=8) :: PSIi,PSIim1,D1n,D1nm1

    !#### Upward recurrence for PSI_n(z) Riccati-Bessel of kind 1
    PSIim1 = ZSIN(z) !i-1=0
    do i=1,n
        call D1nRecu(z,i,nmax,D1n,D1nm1)
        PSIi = PSIim1*(i/z - D1nm1)
        PSIn = PSIi
        PSInm1 = PSIim1
        PSIim1 = PSIi
    enddo

    RETURN
end subroutine PsiRecu

subroutine XiRecu(z,n,nmax,XIn,XInm1)
    implicit none
    !input variables
    INTEGER, INTENT(IN) :: n,nmax
    COMPLEX(kind=8), INTENT(IN) :: z
    !output variables
    COMPLEX(kind=8), INTENT(OUT) :: XIn,XInm1
    !Local variables
    INTEGER :: i
    COMPLEX(kind=8) :: XIi,XIim1,D3n,D3nm1

    !#### Upward recurrence for Xi_n(z) Riccati-Bessel of kind 3
    XIim1 = ZSIN(z) - (0.d0,1.d0)*ZCOS(Z) !i-1=0
    do i=1,n
        call D3nRecu(z,i,nmax,D3n,D3nm1)
        XIi = XIim1*(i/z - D3nm1)
        XIn = XIi
        XInm1 = XIim1
        XIim1 = XIi
    enddo

    RETURN
end subroutine XiRecu

subroutine QRecu(z1,z2,xlm1,xl,n,nmax,Qln)
    implicit none
    !input variables
    INTEGER, INTENT(IN) :: n,nmax
    COMPLEX(kind=8), INTENT(IN) :: z1,z2 !z1=m_(l)x(l-1) and z2=m_(l)x_(l)
    REAL*8 :: xlm1,xl !Size parameter of l-1 and l layers
    !output variables
    COMPLEX(kind=8), INTENT(OUT) :: Qln
    !Local variables
    INTEGER :: i
    REAL*8 :: a1,b1,a2,b2,dn
    COMPLEX(kind=8) :: D1i1,D1i2,D3im11,D3im12,dummy,Qlim1,Qli,dmmy1,dmmy2

    !#### Upward recurrence for Q(l)_n
    a1 = real(z1)
    a1 = -2.*a1
    dmmy1 = COMPLEX(0.d0,a1)
    b1 = imag(z1)
    a2 = real(z2)
    a2 = -2.*a2
    dmmy2 = COMPLEX(0.d0,a2)
    b2 = imag(z2)
    Qlim1 = EXP(dmmy2)-EXP(-2.*b2) !i-1=0
    Qlim1 = (EXP(dmmy1)-EXP(-2.*b1))/Qlim1
    Qlim1 = EXP(-2.*(b2-b1))*Qlim1
    do i=1,n
        dn = i
        call D1nRecu(z1,i,nmax,D1i1,dummy) !D1_i(z1)
        call D1nRecu(z2,i,nmax,D1i2,dummy) !D1_i(z2)
        call D3nRecu(z1,i,nmax,dummy,D3im11) !D3_i(z1)
        call D3nRecu(z2,i,nmax,dummy,D3im12) !D3_i(z2)    
        Qli = (z1*D1i1 + dn)*(dn - z1*D3im11)
        Qli = (z2*D1i2 + dn)*(dn - z2*D3im12)/Qli 
        Qli = Qlim1 * Qli * (xlm1/xl)**2.
        Qlim1 = Qli
    enddo
    Qln = Qli

    RETURN
end subroutine QRecu

real function rads(degrees)
    !Transform degrees into radians
    implicit none
    real,parameter :: PI = 4.0d0*DATAN(1.d0)
    real, intent(in) :: degrees
    
    rads = degrees*PI/180.

end function rads

real function deg(rads)
    !Transform degrees into radians
    implicit none
    real,parameter :: PI = 4.0d0*DATAN(1.d0)
    real, intent(in) :: rads

    deg = rads*180./Pi
    
end function deg

end module MieScattMod