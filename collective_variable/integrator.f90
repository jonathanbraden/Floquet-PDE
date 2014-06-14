#define LAPLACE
#define DVDPHI
! For sine-gordon it's probably faster to simultaneously solve for the breather to avoid evaluating all the stupid atans, etc., but a good idea to check both to ensure agreement
#define D2PHI -cos(4.*atan

#define IFLD 1:nlat
#define IFLDP nlat+1:2*nlat
#define IBG 2*nlat+1:3*nlat
#define IBGP 3*nlat+1:4*nlat
!#define ITIME 4*nlat+1

module integrator
  implicit none

  real(dl) :: k2

contains

! Perform evolution forward in time
!  dt - size of one time step
!  nsteps - number of steps to take (total time evolution is nsteps*dt)
!  fld - values of field to evolve
!  fldp - values of field momentum
!  tglob - current value of time (needed for time-dependent coefficients)
!        - updated as evolution is performed
  subroutine step(dt, nsteps, fld, fldp, tglob)
    real(dl) :: dt
    integer :: nsteps
    

  end subroutine step

!
! The EOM to be solved.  YC is vector of current values, dYC/dt = YP
!
  subroutine derivs(yc, yp)
    real(dl), dimension(1:nvar) :: yc, yp

    real(dl) :: gnorm

    gnorm = 1./dx**2  ! normalization for laplacian

! Currently a very hackish laplacian
!
! Skeleton for spectral laplacian
!  -FFT(fld)
!  -Multiply k^2*fld_k
!  -Inverse FFT into new variable
!
! Requires: Allocating a FT'd array, and allocating an array to store the laplacianed field.  Best to just use in place probably
!

    yp(IFLD) = yc(IFLDP)
    yp(IFLDP) = - (k2 + cos(yc(IBG))) * yc(IFLD) + LAPLACE
    yp(IBG) = yc(IBGP)
    yp(IBGP) = - sin(yc(IBG)) + LAPLACE
!    yp(ITIME) = 1.  ! time variable, only needed if not solving for BG

! D2PHI = V''(phi(t)) 

  end subroutine derivs

!
! The Implicit Runge-Kutta Gauss-Legendre Integrators, the only thing that can be adjusted here is the number of iterative steps
!
  subroutine gl10( y, dt )
    real*8 :: y(nvar)
    real*8 :: dt
    
    integer, parameter :: s = 5
    real*8 :: g(nvar, s)
    
    ! Butcher tableau for 8th order Gauss-Legendre method
    real*8, parameter :: a(s,s) = reshape( (/ &
         0.5923172126404727187856601017997934066D-1, -1.9570364359076037492643214050884060018D-2, &
         1.1254400818642955552716244215090748773D-2, -0.5593793660812184876817721964475928216D-2, &
         1.5881129678659985393652424705934162371D-3,  1.2815100567004528349616684832951382219D-1, &
         1.1965716762484161701032287870890954823D-1, -2.4592114619642200389318251686004016630D-2, &
         1.0318280670683357408953945056355839486D-2, -2.7689943987696030442826307588795957613D-3, &
         1.1377628800422460252874127381536557686D-1,  2.6000465168064151859240589518757397939D-1, &
         1.4222222222222222222222222222222222222D-1, -2.0690316430958284571760137769754882933D-2, &
         4.6871545238699412283907465445931044619D-3,  1.2123243692686414680141465111883827708D-1, &
         2.2899605457899987661169181236146325697D-1,  3.0903655906408664483376269613044846112D-1, &
         1.1965716762484161701032287870890954823D-1, -0.9687563141950739739034827969555140871D-2, &
         1.1687532956022854521776677788936526508D-1,  2.4490812891049541889746347938229502468D-1, &
         2.7319004362580148889172820022935369566D-1,  2.5888469960875927151328897146870315648D-1, &
         0.5923172126404727187856601017997934066D-1 /) , [s,s])
    real, parameter :: b(s) = (/ &
         1.1846344252809454375713202035995868132D-1,  2.3931433524968323402064575741781909646D-1, &
         2.8444444444444444444444444444444444444D-1,  2.3931433524968323402064575741781909646D-1, &
         1.1846344252809454375713202035995868132D-1 /)
    
    integer :: i,k
    
    g = 0.
    do k=1,16
       g = matmul(g,a)
       do i=1,s
          call derivs( y+g(:,i)*dt , g(:,i) )
       enddo
    enddo
    y = y + matmul(g,b)*dt

  end subroutine gl10

end module Integrator
