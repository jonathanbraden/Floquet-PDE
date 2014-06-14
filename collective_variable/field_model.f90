module Model

implicit none

integer, parameter :: dl = kind(1.d0)
real(dl), parameter :: twopi = 6.283185307179586476925867665590_dl
!
! Model parameters, such as lattice size, etc
!
integer, parameter :: nlat = 128 ! Size of the lattice being used'
integer, parameter :: pad = 2  ! padding to implement boundary conditions
! implement this feature in a nicer way
integer, parameter :: imin=1-pad, imax=nlat+pad

!real(dl), dimension(imin:imax) :: fld, fldp
!real(dl), dimension(2) :: rvar
real(dl) :: period
real(dl) :: length
real(dl) :: dx
real(dl) :: x0
real(dl) :: dk

! Model specific parameters
real(dl) :: rmax, epsilon
real(dl), parameter :: width=1./2.**0.5

contains
!
! This is simply to speed up the computations.  We compute the values of the space-dependent coefficients once so we don't have to compute them many times
!
  subroutine set_coeffs(r1, r2)
    real(dl), intent(in) :: r1,r2

! To do, compute r2 rather than specifying as input
    rmax = r1
    epsilon = exp(-2.**1.5*(r1-r2))
    period = 0.5*0.5*twopi*exp(2.*rmax)/6.**0.5  ! compute this via a quadrature integration
    length = 2.*rmax + 40./width
    dx = length / nlat
    x0 = length / 2.
    dk = twopi / length

    print*,"period is ",period," length is ",length
    print*,"dx is ",dx,"x0 is ",x0
  end subroutine set_coeffs

! Adjust this as needed, maybe use phi_max not v
  real(dl) function output_param()
    output_param = period
  end function output_param

! Model independent functions to return dx, period, etc.
  real(dl) function get_period()
    get_period = period
  end function get_period

  real(dl) function get_dx()
    get_dx = dx
  end function get_dx

  real(dl) function k2eff(kreal)
    real(dl) :: kreal

    k2eff = kreal
  end function k2eff

!
! Model specific definition of the effective mass
! Temporary, use a fixed evolution for the background
!
  real(dl) function effective_mass(x,t)
    real(dl), intent(in) :: x,t
    real(dl) :: phitemp
    real(dl), parameter :: winv=1./width

    real(dl), parameter :: rmin=-0.5
    real(dl) :: tmod, xeff, rcur

    tmod = mod(t,period)
    xeff = x-x0
    rcur = rmax + (rmin - rmax)*exp(-(tmod/period-0.5)**2/0.05)

    phitemp = tanh(winv*(xeff+rcur)) - tanh(winv*(xeff-rcur)) - 1
    effective_mass = 3.*phitemp**2 - 1.
  end function effective_mass

  real(dl) function effective_mass_rad(x, rcur)
    real(dl), intent(in) :: x,rcur
    real(dl) :: phitemp, xeff, gamma, width

!    gamma = (1.-rcur(2)**2)**(-0.5)
    width = 1./2.**0.5 ! * gamma
    xeff = x-x0
    phitemp = tanh(width*xeff+rcur) - tanh(width*xeff-rcur) - 1.

    effective_mass_rad = 3.*phitemp**2 - 1.
  end function effective_mass_rad

  real(dl) function r_evolution(t)
    real(dl), intent(in) :: t

    r_evolution = rmax + 0.25*log(cos(0.5*twopi*t/period)**2 + epsilon)
  end function r_evolution

!
! Functions for the background evolution
!
  real(dl) function vprime(r,rdot)
    real(dl), intent(in) :: r, rdot
    real(dl) :: gamma, cth

    gamma = (1.-rdot**2)**(-0.5)
    cth = 1./tanh(2.*r)

    if (abs(r)<1.e-7) then
       vprime = 64.*r/5. - 128.*r**2/5. + 512.*r**3/105. + 512*r**4/21. - 2048.*r**5/175.
! -4096.*r**6/225. + 212992*r**7/17325. + 8192.*r**2/693. - 141197313.*r**9/14189175.
    else
       vprime = 32. - (96.*r+32.)*cth + (96.*r-48.)*cth**2 + (96.*r+48.)*cth**3 - (96.*r*cth**4)
    endif
    vprime = 0.75*vprime
  end function vprime

! I only need this if I want to find r_min numerically
  real(dl) function vdprime(r)
    real(dl), intent(in) :: r
  end function vdprime

  real(dl) function kin_factor(rcur)
    real(dl) :: rcur
    kin_factor = 0.
  end function kin_factor

!
! Gradually increase complexity here
!
  subroutine derivs_bg(yc, yp)
    real, dimension(2), intent(in) :: yc
    real, dimension(2), intent(out) :: yp

    real(dl) :: gamma
    gamma = (1.-yc(2)**2)**0.5
    yp(1) = yc(2)
    yp(2) = -vprime(yc(1),yc(2)) !/ ( gamma**3 ) !* 2. * tension )
  end subroutine derivs_bg

end module Model
