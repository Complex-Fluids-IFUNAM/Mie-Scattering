module inputdata

    implicit none
    !##############################################################################!
    !#######                         INPUT DATA MODULE                      #######!
    !####### MIT License                                                    #######!        
    !####### Copyright (c) 2024 Christian Balderas, Rolando Castillo &      #######!
    !####### Complex-Fluids-IFUNAM                                          #######!
    !##############################################################################!      

    real :: WAVEL,REFMED,k0,Rp
    complex :: np
    integer :: NANG
    logical :: HOMOSPHE
    !REAL,ALLOCATABLE,DIMENSION(:) :: xl
    !complex,ALLOCATABLE,DIMENSION(:) :: ml

    CONTAINS
        subroutine inputs(WAVEL,REFMED,k0,HOMOSPHE,NANG)
            implicit none
            REAL,INTENT(OUT) :: WAVEL,REFMED,k0
            logical, INTENT(OUT) :: HOMOSPHE
            integer, INTENT(OUT) :: NANG

            REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
            
            NAMELIST / PARPHYS / WAVEL,REFMED,HOMOSPHE,NANG

            OPEN(UNIT=1,FILE='parameters.nml',DELIM='APOSTROPHE')
            READ(UNIT=1,nml=PARPHYS)
            CLOSE(UNIT=1)

            k0 = 2*PI*REFMED/WAVEL          !wavenumber k0

            RETURN
        end subroutine inputs

        !subroutine AllocLyrsArrays(ly)
        !    implicit none
        !    INTEGER, INTENT(IN) :: ly
        
            !Allocate memory for size and RI arrays
        !    if (.NOT. ALLOCATED(xl)) then
        !        ALLOCATE(xl(ly),ml(ly))
        !    ELSE
        !        DEALLOCATE(xl,ml)
        !        ALLOCATE(xl(ly),ml(ly))
        !    endif
    
        !end subroutine AllocLyrsArrays

        subroutine ReadRIprofile(Aly,xl,ml)
            implicit none
            !input variables
            INTEGER, INTENT(IN) :: Aly
            !output variables 
            REAL, INTENT(OUT) :: xl(Aly)
            complex, INTENT(OUT) :: ml(Aly)
            !Local variables
            integer :: i
            
            !#### READ RI profile table ######!
            !#### WATCH UP!!!!!!!!!!!
            !#### The array must have columns [Radius,RI(Radius)] with Radius in nm.
            !#### Numerical format -------->  [ Real ,  (Re,Im) ] 
            !#### If the table is hasn't this formar You must take action to guarantee this.
            OPEN(UNIT=1,FILE='RIprofile.dat',STATUS='OLD',ACTION='READ')

            !#### Reads RI profile data ####!
            do i=1,Aly
                READ(1,*) xl(i),ml(i)
            enddo
            CLOSE(UNIT=1)
        
            RETURN
        end subroutine ReadRIprofile

end module inputdata