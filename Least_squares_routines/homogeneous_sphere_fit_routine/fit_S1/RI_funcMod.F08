module RI_funcMod
    implicit none
    !#########################################################
    !#### This module contains functions and subroutines  ####
    !#### for estimating the refractive index of a sphere ####
    !#### with a particular behavior.                     ####
    !#### C. Balderas, August 2022.                       ####
    !#### Complex Fluids, UNAM-MX                         ####
    !#########################################################
    INTEGER :: Aly,ly
    REAL,ALLOCATABLE,DIMENSION(:) :: xl
    COMPLEX,ALLOCATABLE,DIMENSION(:) :: ml

    CONTAINS
    
    subroutine AllocLyrsArrays(ly)
        implicit none
        INTEGER, INTENT(IN) :: ly
    
        !Allocate memory for size and RI arrays
        if (.NOT. ALLOCATED(xl)) then
            ALLOCATE(xl(ly),ml(ly))
        ELSE
            DEALLOCATE(xl,ml)
            ALLOCATE(xl(ly),ml(ly))
        endif

    end subroutine AllocLyrsArrays

    subroutine CSGaussProf(Rcos,Rs,REFMED,RImax,sigm,xl,ml,ly)
        !use inputdata
        implicit none
        !## This subroutine calculates a gaussian refractive index profile 
        !Input and output variables
        REAL, INTENT(in) :: REFMED,Rs,Rcos,sigm !Rcos = Core of the shell Radius
        COMPLEX, INTENT(in) :: RImax
        REAL,INTENT(inout) :: xl(Aly)
        COMPLEX,INTENT(inout) :: ml(Aly)
        INTEGER,INTENT(out) :: ly      
        !Local variables
        INTEGER :: i,j,l
        REAL :: dummy,rm,dR
        COMPLEX :: RI
        
        !FINDING CUT RADIUS, MUST BE <= HYDRODINAMIC RADIUS R_h
        dR = 0.1
        RM = Rcos
        RI = RImax
        i = 1
        do while (abs(RI-REFMED) .GE. 0.0001 .AND. RM .LE. Rs)
            RM = Rcos + real(i)*DR
            RI = REFMED + (RImax-REFMED)*halfgauss(RM,Rcos,sigm)
            i = i+1
        enddo
        i = i-1
        if (i .LT. Aly) then
            ly = i
        else
            ly = Aly
        endif
        DR = (RM-Rcos)/ly
        !CALC SIZE AND RI OF i-layer - i=1 is reserved for solid core
        do i=2,ly
            xl(i) = Rcos + real(i-1)*DR
        enddo
        ml(2) = RImax
        do i=3,ly
            rm = (xl(i)+xl(i-1))*0.5
            ml(i) = REFMED + (RImax-REFMED)*halfgauss(RM,Rcos,sigm)
        enddo

        RETURN
    end subroutine CSGaussProf

    subroutine GaussProf(Rcore,Rcut,REFMED,RImax,sigm,xl,ml,ly)
        implicit none
        !## This subroutine calculates a gaussian refractive index profile 
        !Input and output variables
        REAL, INTENT(in) :: REFMED,Rcut,Rcore,sigm
        COMPLEX, INTENT(in) :: RImax
        REAL,INTENT(inout) :: xl(Aly)
        COMPLEX,INTENT(inout) :: ml(Aly)
        INTEGER,INTENT(out) :: ly      
        !Local variables
        INTEGER :: i,j,l
        REAL :: dummy,rm,dR
        COMPLEX :: RI
        
        !FINDING CUT RADIUS, MUST BE <= HYDRODINAMIC RADIUS R_h
        dR = 0.1
        RM = Rcore
        RI = RImax
        i = 1
        do while (abs(RI-REFMED) .GE. 0.0001 .AND. RM .LE. Rcut)
            RM = Rcore + real(i)*DR
            RI = REFMED + (RImax-REFMED)*halfgauss(RM,Rcore,sigm)
            i = i+1
        enddo
        i = i-1
        if (i .LT. Aly) then
            ly = i
        else
            ly = Aly
        endif
        DR = (RM-Rcore)/ly
        !CALC SIZE AND RI OF i-layer
        do i=1,ly
            xl(i) = Rcore + real(i-1)*DR
        enddo
        ml(1) = RImax
        do i=2,ly
            rm = (xl(i)+xl(i-1))*0.5
            ml(i) = REFMED + (RImax-REFMED)*halfgauss(RM,Rcore,sigm)
        enddo

        !OPEN(UNIT=1,FILE='ri.dat',STATUS='UNKNOWN',ACTION='WRITE')
        !do i=1,ly
        !    write(1,*) xl(i),real(ml(i))
        !    write(*,*) xl(i),real(ml(i))
        !enddo
        !close(unit=1)
        !stop

        RETURN
    end subroutine GaussProf

    real function halfgauss(x,x0,sigm)
        implicit none
        real, INTENT(in) :: x,x0,sigm
        real :: dummy

        dummy = (x-x0)/sigm
        halfgauss = EXP(-0.5 * dummy**2.)

    end function halfgauss

    real function sigmgauss(Rs,Rc,nc,ns,n0)
        implicit none
        real, INTENT(in) :: Rs,Rc,nc,ns,n0
        real :: dummy

        dummy = (nc - n0)/(ns - n0)
        sigmgauss = (Rs-Rc)/sqrt(2.*log(dummy))

    end function sigmgauss

end module RI_funcMod
