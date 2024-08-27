module inputdata

    implicit none
    real :: WAVEL,REFMED,k0,HRAD,Rps,stdgps
    complex :: nps
    integer :: NANG
    REAL,ALLOCATABLE,DIMENSION(:,:) :: data

    CONTAINS
        subroutine inputs(WAVEL,REFMED,k0,HRAD)!,Rps,stdgps,nps)
            implicit none
            REAL,INTENT(OUT)::WAVEL,REFMED,k0,HRAD!,Rps,stdgps
            !complex,INTENT(out)::nps
            REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
            
            NAMELIST / PARPHYS / WAVEL,REFMED,HRAD
            OPEN(UNIT=1,FILE='inputs.nml',DELIM='APOSTROPHE')
            READ(UNIT=1,nml=PARPHYS)
            CLOSE(UNIT=1)
            !WAVEL = 632.8                   !Wavelenght (nm) 
            !REFMED = 1.332                  !refractive index of medium
            k0 = 2*PI*REFMED/WAVEL          !wavenumber k0
            !Rps = 0.                      !PS particle radius (nm)
            !nps = (0.,0.)                 !PS refractive index
            !stdgps = 0.                   !PS particle polydispersity
            !HRAD = 180.                     !Hydrodinamic radius (nm)

            RETURN
        end subroutine inputs

        subroutine ReadExpData(NANG,data)
            implicit none
            REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
            !output variables 
            integer, intent(out) :: NANG 
            real,allocatable,dimension(:,:),intent(out) :: data
            integer :: i,status
            real :: dummy
            
            !#### READ EXPERIMENTAL DATA ######!
            !#### WATCH UP!!!!!!!!!!!
            !#### The array must have columns [theta,I(theta),ErrBars] with theta in degrees. 
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
                !### ONLY USE WHEN I=I(q), theta=2*ASIN(q*lmda/2*k0) ###! 
                !data(1,i) = 180.*(2.*ASIN(data(1,i)/(2*k0)))/PI 
            enddo
            CLOSE(UNIT=1)
        
            RETURN
        end subroutine ReadExpData

end module inputdata