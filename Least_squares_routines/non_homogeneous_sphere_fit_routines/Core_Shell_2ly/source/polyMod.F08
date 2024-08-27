module polyMod
implicit none
!##############################################################################!
!#######                POLYDISPERSITY INTEGRALS MODULE                 #######!
!####### This module contain the numerical method to calculate          #######!
!####### size polydispersity integrals. Gauss-Legendre quadrature is    #######!
!####### used. The main algorithm was taken from:                       #######!
!####### Teukolsky, S. A., & Vetterling, W. T. (1992).                  #######!
!####### Numerical Methods for Fortran 77: The art of scientific        #######!
!####### computing.                                                     #######!
!####### MIT License                                                    #######!        
!####### Copyright (c) 2024 Christian Balderas, Rolando Castillo &      #######!
!####### Complex-Fluids-IFUNAM                                          #######!
!##############################################################################!      
REAL :: QSCApoly
REAL,ALLOCATABLE,DIMENSION(:) :: NS1poly,NS2poly,Fpoly

contains
    subroutine AllocPolyArray(NANG)
        implicit none
        INTEGER,intent(IN) :: NANG
    
        if (.not. allocated(NS1poly)) then
            ALLOCATE(NS1poly(NANG),NS2poly(NANG),Fpoly(NANG))
        else
            DEALLOCATE(NS1poly,NS2poly,Fpoly)
            ALLOCATE(NS1poly(NANG),NS2poly(NANG),Fpoly(NANG))
        endif
        
        RETURN
    end subroutine AllocPolyArray
!
    subroutine polyscatt(stdgfit,QSCApoly,Fpoly,NS1poly,NS2poly)
    !##############################################################################
    !### POLYDISPERSITY SUBROUTINE FOR A SINGLE SPHERE                          ###
    !##############################################################################
        USE inputdata
        USE MieScattMod
        implicit none
        
        !## INOUT VARIABLES ##!
        INTEGER, PARAMETER :: N=8 !LEGENDRE POLYNOMIAL OF N-th order
        REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
        REAL, INTENT(in) :: stdgfit
        REAL,INTENT(OUT) :: QSCApoly,Fpoly(NANG),NS1poly(NANG),NS2poly(NANG)
        !## LOCAL VARIABLES ##!
        INTEGER i,j,k
        REAL :: a,b,R,dR,xm,xr,w(n),t(n),&
                QSCAm
        !##########################################################################
        !Parametrization:                                                          
        !R(t) = xm + xr*t
        !dR=xr*dt
        !Gauss quadrature: int[f(x)dx]_[a,b]=int[f(xm + xr*t)xr*dt]_[-1,1]
        !                                   =xr*Sum[w_i * f(xm + xr*t_i)]
        !##########################################################################
        a = X*(stdgfit**(-sqrt(18.))) !Lower integration limit
        b = X*(stdgfit**(sqrt(18.)))  !Upper integration limit
        call gauleg(-1.,1.,n,t,w) !Roots of legendre polynom with order n
        xm = 0.5*(b+a)  !Central point in [a,b] 
        xr = 0.5*(b-a)  !slope in parametrization
        QSCApoly = 0.   !quadrature initialization
        Fpoly(:) = 0.
        NS1poly(:) = 0.
        NS2poly(:) = 0. 
        do j=1,n !CROSS SECTION WITH POLYDISPERSITY
            dR=xr*t(j)
            R=xm+dR
            CALL BHMIE(R,REFREL,NS1,NS2,FF,QSCAm)
            QSCApoly = QSCApoly + w(j)*QSCAm*lgnormfun(R,X,stdgfit)
            Fpoly(:) = Fpoly(:) + w(j)*FF(:)*lgnormfun(R,X,stdgfit)
            NS1poly(:) = NS1poly(:) + w(j)*NS1(:)*lgnormfun(R,X,stdgfit)
            NS2poly(:) = NS1poly(:) + w(j)*NS2(:)*lgnormfun(R,X,stdgfit)
        enddo
        QSCApoly = xr*QSCApoly*(exp(-0.5*(log(stdgfit))**2.))
        Fpoly(:) = xr*Fpoly(:)*(exp(-0.5*(log(stdgfit))**2.))    
        NS1poly(:) = xr*NS1poly(:)*(exp(-0.5*(log(stdgfit))**2.))
        NS2poly(:) = xr*NS1poly(:)*(exp(-0.5*(log(stdgfit))**2.))
        RETURN
    end subroutine polyscatt

    subroutine TwoLys_PolyRs(Rc,Rs,stdgfit,QSCApoly,Fpoly,NS1poly,NS2poly)
    !##############################################################################
    !### POLYDISPERSITY SUBROUTINE FOR MULTI-LAYERED SPHERE WITH MONODISPERSE   ###
    !### CORE RADIUS. THE REFRACTIVE INDEX HAS A GAUSSIAN PROFILE WITH r        ###
    !##############################################################################
    !Parametrization:                                                           ###
    !R(t) = xm + xr*t                                                           ###
    !dR=xr*dt                                                                   ###
    !Gauss quadrature: int[f(x)dx]_[a,b]=int[f(xm + xr*t)xr*dt]_[-1,1]          ###
    !                                   =xr*Sum[w_i * f(xm + xr*t_i)]           ###
    !##############################################################################
    !### BY C. BALDERAS.   2022                                                 ###
    !##############################################################################
        USE inputdata
        USE MieScattMod
        USE RI_funcMod
        USE omp_lib
        implicit none
        !## INOUT VARIABLES ##!
        INTEGER, PARAMETER :: N=8 !LEGENDRE POLYNOMIAL OF N-th order
        REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
        REAL, INTENT(in) :: stdgfit,Rc,Rs
        !COMPLEX,INTENT(in) :: RIfit
        REAL,INTENT(OUT) :: QSCApoly,Fpoly(NANG),NS1poly(NANG),NS2poly(NANG)
        !## LOCAL VARIABLES ##!
        INTEGER i,j,TID
        REAL :: a,b,R,dR,xm,xr,w(n),t(n),&
                QEXTm,QSCAm,QBACKm,CSCAm,&
                Xs,Xc,sigm,ratio,ti,tf
        !ratio = Rs/Rcore
        Xs = Rs*k0
        a = Xs*(stdgfit**(-sqrt(18.))) !Lower integration limit
        b = Xs*(stdgfit**(sqrt(18.)))  !Upper integration limit
        call gauleg(-1.,1.,n,t,w) !Roots of legendre polynom with order n
        xm = 0.5*(b+a)  !Central point in [a,b] 
        xr = 0.5*(b-a)  !slope in parametrization
        QSCApoly = 0.   !quadrature initialization
        Fpoly(:) = 0.
        NS1poly(:) = 0.
        NS2poly(:) = 0.
        !######## THIS SECTION IS FOR A PARALLEL THREADS WORK USING OPEN_MP ########################!
        !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,dR,R,sigm,xl,ml,ly,QEXTm,QSCAm,QBACKm,FF,NS1,NS2) 
        !TID = OMP_GET_THREAD_NUM()
        !!$OMP DO PRIVATE(j,i)
        do j=1,n !CROSS SECTION WITH POLYDISPERSITY
            dR=xr*t(j)
            R=xm+dR
            xl(2) = R
            CALL MieMLy(ly,ml,xl,QEXTm,QSCAm,QBACKm,FF,NS1,NS2)
            !!$OMP CRITICAL (poly_arrays) !Output actualization
            QSCApoly = QSCApoly + w(j)*QSCAm*lgnormfun(R,Xs,stdgfit)
            Fpoly(:) = Fpoly(:) + w(j)*FF(:)*lgnormfun(R,Xs,stdgfit)
            NS1poly(:) = NS1poly(:) + w(j)*NS1(:)*lgnormfun(R,Xs,stdgfit)
            NS2poly(:) = NS2poly(:) + w(j)*NS2(:)*lgnormfun(R,Xs,stdgfit)
            !!$OMP END CRITICAL (poly_arrays)
        enddo
        !!$OMP END DO
        !!$OMP END PARALLEL
        QSCApoly = xr*QSCApoly*(exp(-0.5*(log(stdgfit))**2.))
        Fpoly(:) = xr*Fpoly(:)*(exp(-0.5*(log(stdgfit))**2.))
        NS1poly(:) = xr*NS1poly(:)*(exp(-0.5*(log(stdgfit))**2.))
        NS2poly(:) = xr*NS1poly(:)*(exp(-0.5*(log(stdgfit))**2.))
        
        RETURN
    end subroutine TwoLys_PolyRs

    subroutine TwoLys_PolyRc(Rc,Rs,stdgfit,QSCApoly,Fpoly,NS1poly,NS2poly)
    !##############################################################################
    !### POLYDISPERSITY SUBROUTINE FOR MULTI-LAYERED SPHERE WITH MONODISPERSE   ###
    !### CORE RADIUS. THE REFRACTIVE INDEX HAS A GAUSSIAN PROFILE WITH r        ###
    !##############################################################################
    !Parametrization:                                                           ###
    !R(t) = xm + xr*t                                                           ###
    !dR=xr*dt                                                                   ###
    !Gauss quadrature: int[f(x)dx]_[a,b]=int[f(xm + xr*t)xr*dt]_[-1,1]          ###
    !                                   =xr*Sum[w_i * f(xm + xr*t_i)]           ###
    !##############################################################################
    !### BY C. BALDERAS.   2022                                                 ###
    !##############################################################################
        USE inputdata
        USE MieScattMod
        USE RI_funcMod
        USE omp_lib
        implicit none
        !## INOUT VARIABLES ##!
        INTEGER, PARAMETER :: N=8 !LEGENDRE POLYNOMIAL OF N-th order
        REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
        REAL, INTENT(in) :: stdgfit,Rc,Rs
        !COMPLEX,INTENT(in) :: RIfit
        REAL,INTENT(OUT) :: QSCApoly,Fpoly(NANG),NS1poly(NANG),NS2poly(NANG)
        !## LOCAL VARIABLES ##!
        INTEGER i,j,TID
        REAL :: a,b,R,dR,xm,xr,w(n),t(n),&
                QEXTm,QSCAm,QBACKm,CSCAm,&
                Xs,Xc,sigm,ratio,ti,tf
        !ratio = Rs/Rcore
        Xc = Rc*k0
        a = Xc*(stdgfit**(-sqrt(18.))) !Lower integration limit
        b = Xc*(stdgfit**(sqrt(18.)))  !Upper integration limit
        call gauleg(-1.,1.,n,t,w) !Roots of legendre polynom with order n
        xm = 0.5*(b+a)  !Central point in [a,b] 
        xr = 0.5*(b-a)  !slope in parametrization
        QSCApoly = 0.   !quadrature initialization
        Fpoly(:) = 0.
        NS1poly(:) = 0.
        NS2poly(:) = 0.
        do j=1,n !CROSS SECTION WITH POLYDISPERSITY
            dR=xr*t(j)
            R=xm+dR
            xl(1) = R
            CALL MieMLy(ly,ml,xl,QEXTm,QSCAm,QBACKm,FF,NS1,NS2)
            QSCApoly = QSCApoly + w(j)*QSCAm*lgnormfun(R,Xs,stdgfit)
            Fpoly(:) = Fpoly(:) + w(j)*FF(:)*lgnormfun(R,Xs,stdgfit)
            NS1poly(:) = NS1poly(:) + w(j)*NS1(:)*lgnormfun(R,Xs,stdgfit)
            NS2poly(:) = NS2poly(:) + w(j)*NS2(:)*lgnormfun(R,Xs,stdgfit)
        enddo
        QSCApoly = xr*QSCApoly*(exp(-0.5*(log(stdgfit))**2.))
        Fpoly(:) = xr*Fpoly(:)*(exp(-0.5*(log(stdgfit))**2.))
        NS1poly(:) = xr*NS1poly(:)*(exp(-0.5*(log(stdgfit))**2.))
        NS2poly(:) = xr*NS1poly(:)*(exp(-0.5*(log(stdgfit))**2.))
        
        RETURN
    end subroutine TwoLys_PolyRc

    subroutine TwoLys_PolyRcRs(Rc,Rs,stdgRc,stdgRs,nc,ns,QSCApoly,Fpoly,NS1poly,NS2poly)
    !##############################################################################
    !### POLYDISPERSITY SUBROUTINE FOR MULTI-LAYERED SPHERE WITH POLYDISPERSITY ###
    !### IN BOTH CORE AND SHELL RADIUS.                                         ###
    !### THE REFRACTIVE INDEX HAS A GAUSSIAN PROFILE WITH r                     ###
    !##############################################################################
    !Parametrization:                                                           ###
    !R(t) = xm + xr*t                                                           ###
    !dR=xr*dt                                                                   ###
    !Gauss quadrature: int[f(x)dx]_[a,b]=int[f(xm + xr*t)xr*dt]_[-1,1]          ###
    !                                   =xr*Sum[w_i * f(xm + xr*t_i)]           ###
    !##############################################################################
    !### BY C. BALDERAS.   2022                                                 ###
    !##############################################################################
        USE inputdata
        USE MieScattMod
        USE RI_funcMod
        USE omp_lib
        implicit none
        !## INOUT VARIABLES ##!
        INTEGER, PARAMETER :: N=8 !LEGENDRE POLYNOMIAL OF N-th order
        REAL, PARAMETER :: PI = 4.0d0*DATAN(1.d0) !Pi=3.1415
        REAL, INTENT(in) :: Rc,Rs,stdgRc,stdgRs
        COMPLEX,INTENT(in) :: nc,ns
        REAL,INTENT(OUT) :: QSCApoly,Fpoly(NANG),NS1poly(NANG),NS2poly(NANG)
        !## LOCAL VARIABLES ##!
        INTEGER i,j,TID
        REAL :: ax,bx,xm,xr,p,dp,Xc,&
                ay,by,ym,yr,q,dq,Xs,w(n),t(n),&
                QSCAsi,Fsi(NANG),NS1si(NANG),NS2si(NANG),&
                QEXTm,QSCAm,QBACKm,CSCAm,&
                sigm,tf,ti
        !------------ Set abssisas and parametrization ------------!
        !------------ x-axis
        Xc = Rc*k0
        ax = Xc*(stdgRc**(-sqrt(18.)))  !Lower integration limit
        bx = Xc*(stdgRc**(sqrt(18.)))   !Upper integration limit
        xm = 0.5*(bx+ax)  !Central point in [ax,bx] 
        xr = 0.5*(bx-ax)  !slope in parametrization
        !------------ y-axis
        Xs = Rs*k0
        ay = Xs*(stdgRs**(-sqrt(18.)))  !Lower integration limit
        by = Xs*(stdgRs**(sqrt(18.)))   !Upper integration limit
        ym = 0.5*(by+ay)  !Central point in [ay,by] 
        yr = 0.5*(by-ay)  !slope in parametrization
        !----------------------------------------------------------!
        !write(*,*) Rps,ax/k0,bx/k0,Rc,ay/k0,by/k0,"Rs",Rs
        if (bx .GE. ay) then
            write(*,*) "The R shell cannot be less than Rcore"
            write(*,*) "stdgRc =",stdgRc,"stdgRs=",stdgRs
            write(*,*) "Rc_max = ",bx/k0,"Rs_min=",ay/k0 
            stop
        endif
        call gauleg(-1.,1.,n,t,w) !Roots of legendre polynom with order n
        !quadrature initialization
        QSCApoly = 0.
        Fpoly = 0.
        NS1poly = 0.
        NS2poly = 0.
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID,dq,q,dp,p,xl,ml,Aly,QEXTm,QSCAm,QBACKm,FF,NS1,NS2,QSCAsi,Fsi,NS1si,NS2si) 
        TID = OMP_GET_THREAD_NUM()
        Aly = 2
        xl(1) = Rc
        xl(2) = Rs
        xl(:) = k0*xl(:)
        ml(1) = nc
        ml(2) = ns
        ml(:) = ml(:)/REFMED
        !$OMP DO PRIVATE(i,j)
        do j=1,n !CROSS SECTION WITH POLYDISPERSITY
            dq=yr*t(j)
            q=ym+dq
            QSCAsi = 0.
            Fsi(:) = 0.
            NS1si(:) = 0.
            NS2si(:) = 0.
            do i=1,n
                dp=xr*t(i)
                p=xm+dp
                xl(1) = p
                xl(2) = q
                CALL MieMLy(Aly,ml,xl,QEXTm,QSCAm,QBACKm,FF,NS1,NS2)
                !write(*,*) "thread",TID,"Qsca=",QSCAm,"Aly=",Aly
                QSCAsi = QSCAsi + w(i)*QSCAm*lgnormfun(p,Xc,stdgRc)*lgnormfun(q,Xs,stdgRs)
                Fsi(:) = Fsi(:) + w(i)*FF(:)*lgnormfun(p,Xc,stdgRc)*lgnormfun(q,Xs,stdgRs)
                NS1si(:) = NS1si(:) + w(i)*NS1(:)*lgnormfun(p,Xc,stdgRc)*lgnormfun(q,Xs,stdgRs)
                NS2si(:) = NS2si(:) + w(i)*NS2(:)*lgnormfun(p,Xc,stdgRc)*lgnormfun(q,Xs,stdgRs)
            enddo
            !WRITE(*,'(A,I2,A,F6.2,A,ES14.6)') 'Hello from thread ',TID!,', I take X=',X,'and the Qsca is',QSCApoly
            !$OMP CRITICAL (poly_arrays)
            QSCApoly = QSCApoly + w(j)*QSCAsi
            !write(*,*) "thread",TID,"Qsca=",QSCApoly
            Fpoly(:) = Fpoly(:) + w(j)*Fsi(:)
            NS1poly(:) = NS1poly(:) + w(j)*NS1si(:)
            NS2poly(:) = NS2poly(:) + w(j)*NS2si(:)
            !$OMP END CRITICAL (poly_arrays)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        !tf = OMP_GET_WTIME()
        !write(*,*) "Work in parallel time = ",tf-ti
        QSCApoly = xr*yr*QSCApoly*(exp(-0.5*(log(stdgRc))**2.))*(exp(-0.5*(log(stdgRs))**2.))
        Fpoly(:) = xr*yr*Fpoly(:)*(exp(-0.5*(log(stdgRc))**2.))*(exp(-0.5*(log(stdgRs))**2.))
        NS1poly(:) = xr*yr*NS1poly(:)*(exp(-0.5*(log(stdgRc))**2.))*(exp(-0.5*(log(stdgRs))**2.))
        NS2poly(:) = xr*yr*NS1poly(:)*(exp(-0.5*(log(stdgRc))**2.))*(exp(-0.5*(log(stdgRs))**2.))
        RETURN
    end subroutine TwoLys_PolyRcRs

    real function lgnormfun(x,x0,stdg)
        implicit none
        !############################################################
        !Log normal function:                                       #
        !                1                                          #
        ! f(x)=_____________________ exp{-1/2[ln(x/x0)/ln(stdg)]^2} #
        !       x (2Pi)^0.5 ln(stdg)                                #
        !                                                           #
        !############################################################
        REAL, INTENT(IN) :: x,x0,stdg
        !REAL, INTENT(OUT) :: lgnorm
        REAL :: d1,d2,d3,PI
        Pi=4.0d0*DATAN(1.d0)
        !Evaluate expression
        d1 = (LOG(x/x0)/LOG(stdg))**2.
        d2 = EXP(-0.5*d1)
        d3 = 1/(SQRT(2*Pi)*LOG(stdg))
        lgnormfun = d2*d3/x
        
    end function lgnormfun
    
    subroutine gauleg(x1,x2,n,x,w)
        implicit none
        !INOUT VARIABLES
        REAL,INTENT(IN) :: x1,x2
        INTEGER,INTENT(IN) :: n
        REAL,INTENT(OUT) :: x(n),w(n)
        !INTERNAL VARIABLES
        REAL*8 :: EPS,p1,p2,p3,pp,xl,xm,z,z1,pi,absc,weight
        INTEGER :: i,j,k,m
        PI = 4.0d0*DATAN(1.d0)
        EPS = 3.d-16
        !################################
        m = Int((n+1)/2) !The roots are symmetric, so we only have to find half of them
        k = n
        xm = 0.5d0*(real(x2,8)+real(x1,8)) !Parametrization
        xl = 0.5d0*(real(x2,8)-real(x1,8))
        do i=1,m
            z = dcos(pi*(real(i,8)-0.25d0)/(real(n,8)+0.5d0))
            !write(*,*) pi*(real(i,8)-0.25d0)/(real(n,8)+0.5d0)  !debugger
            !pause (1)
            do while (dabs(z-z1) .ge. EPS) !Finding root with Newton's method
                p1 = 1.d0
                p2 = 0.d0
                do j=1,n
                    !write(*,*) "inner loop", j
                    p3=p2
                    p2=p1
                    p1=((2.d0*real(j,8)-1.d0)*z*p2-(real(j,8)-1.d0)*p3)/real(j,8)
                    !p1 is now the desired Legendre Polynomial
                enddo
                !Recurrence relations for the Legendre polynomials
                pp=real(k,8)*(z*p1-p2)/(z*z-1.d0)
                z1 = z
                z = z1 - p1/pp 
            enddo
            x(i) = xm - xl*z        !Parametrization
            x(n+1-i) = xm + xl*z    !Symmetric counterpart
            w(i) = 2.d0*xl/((1.d0 - z**2.d0)*pp**2.d0)        !Compute the weight
            w(n+1-i) = w(i)                                   !and its symmetric counterpart
            !write(*,*) "debugger,i=",i,"x(i)=",x(i),"w(i)=",w(i),"x(k+1-i)=",x(k+1-i),"w(k+1-i)=",w(k+1-i)
        enddo
    end subroutine gauleg
    
    subroutine qgauss2d(integral)
        !##############################################################################!
        !#### THIS SUBROUTINE IS A TEST FOR THE GAUSSIAN-LEGENDRE IN 2D ALGORITHM #####!
        !#### THE PROVE FUNCTION IS: x^4*y^2 + 2*y^2*x^2 - y*x^2 + 2              #####!
        !#### THE INTEGRAL F(x,y)dxdy over x=[-1,1], y=[-1,1] is 9.15556 (exact)  #####!
        !##############################################################################!
        implicit none
        INTEGER, PARAMETER :: n = 5 !Legendre polynom with order n
        real, intent(out) :: integral
        integer :: i,j
        real :: ax,bx,xm,xr,p,dp,&
                ay,by,ym,yr,q,dq,&
                w(n),t(n),sumi

        !----------- Set the absissas -------------!
        !#### x-axis
        ax = -1.
        bx = 1.
        xm = 0.5*(bx+ax)  !Central point in [a,b] 
        xr = 0.5*(bx-ax)  !slope in parametrization
        !#### y-axis
        ay = -1.
        by = 1.
        ym = 0.5*(by+ay)
        yr = 0.5*(by-ay)
        !------------------------------------------!
        call gauleg(-1.,1.,n,t,w) !Roots of legendre polynom with order n
        integral = 0. !Quadrature initialization
        do j=1,n
            dq=yr*t(j)
            q=ym+dq
            sumi = 0. 
            do i=1,n
                dp=xr*t(i)
                p=xm+dp
                sumi = sumi + w(i)*provf(p,q)
            enddo
            integral = integral + w(j)*sumi
        enddo

        return
    end subroutine qgauss2d

    real function provf(x,y)
        implicit none
        real,intent(in) :: x,y

        provf = (x**4.)*(y**2.) + 2*(y**2.)*(x**2.) - y*(x**2.) + 2.

    end function provf

end module polyMod


