module inputdata

    implicit none
    !##############################################################################!
    !#######                         INPUT DATA MODULE                      #######!
    !####### ONLY FOR RI PROFILE MODELED AS: CORE-RI PLATEAU-GAUSSIAN DECAY #######!
    !####### This module contain two subroutines. The first one reads the   #######!
    !####### experimental parameters like wavelenght, core's RI and         #######!
    !####### polydispersity.                                                #######!
    !####### MIT License                                                    #######!        
    !####### Copyright (c) 2024 Christian Balderas, Rolando Castillo &      #######!
    !####### Complex-Fluids-IFUNAM                                          #######!
    !##############################################################################!      

    real :: WAVEL,REFMED,k0,HRAD,CHI2MIN
    integer :: NANG
    REAL,ALLOCATABLE,DIMENSION(:,:) :: data

    CONTAINS
        subroutine inputs(WAVEL,REFMED,k0,HRAD,CHI2MIN)
            implicit none
            REAL,INTENT(OUT)::WAVEL,REFMED,k0,HRAD,CHI2MIN
            REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
            
            NAMELIST / PARPHYS / WAVEL,REFMED,HRAD,CHI2MIN
            OPEN(UNIT=1,FILE='parameters.nml',DELIM='APOSTROPHE')
            READ(UNIT=1,nml=PARPHYS)
            CLOSE(UNIT=1)

            k0 = 2*PI*REFMED/WAVEL          !wavenumber k0

            RETURN
        end subroutine inputs

        subroutine ReadExpData(NANG,data)
            implicit none
            !output variables 
            integer, intent(out) :: NANG 
            real,allocatable,dimension(:,:),intent(out) :: data
            integer :: i,status
            real :: dummy
            
            !#### READ EXPERIMENTAL DATA ######!
            !#### WATCH UP!!!!!!!!!!!
            !#### The array must have columns [theta,I(theta)] with theta in degrees. 
            !#### If the data is not a function of theta or is not in degrees, 
            !#### You must take action to guarantee this.
            OPEN(UNIT=1,FILE='intensity.dat',STATUS='OLD',ACTION='READ')
            !## what's the experimental data array size?
            i = 0
            status = 0
            do while (status == 0)
                READ(1,*,IOSTAT=status) dummy, dummy, dummy !
                i = i + 1
            enddo
            NANG = i - 1
            !#### Allocate memory for data array ####!
            if (.NOT. allocated(data)) then
                allocate(data(3,NANG))
            else
                deallocate(data)
                allocate(data(3,NANG))
            endif
            !#### Write experimental data on data array ####!
            REWIND(UNIT=1)
            do i=1,NANG
                READ(1,*) data(1,i), data(2,i), data(3,i)
            enddo
            CLOSE(UNIT=1)
        
            RETURN
        end subroutine ReadExpData

end module inputdata