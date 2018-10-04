module Model
  use constants
  implicit none
!
! Model parameters, such as lattice size, etc
!
integer, parameter :: nlat = 32  ! Size of the lattice being used'
integer, parameter :: pad = 2  ! padding to implement boundary conditions
! implement this feature in a nicer way
integer, parameter :: imin=1-pad, imax=nlat+pad

! make these allocatable
real(dl), dimension(imin:imax) :: fld, fldp, fldbg, fldpbg
real(dl) :: period
real(dl) :: length
real(dl) :: dx
real(dl) :: x0
real(dl) :: dk

! Model specific parameters
real(dl) :: epsilon
real(dl) :: proamp
real(dl) :: omega
real(dl) :: width
real(dl) :: geff
real(dl),parameter :: rat=1.e-6  ! value of background breather at boundaries

contains
!
! This is simply to speed up the computations.  We compute the values of the space-dependent coefficients once so we don't have to compute them many times
!
  subroutine set_coeffs(epsset)
    real(dl) :: epsset

    epsilon = epsset

    proamp = 2.*epsilon / sqrt(3.)
! Set this to 0 to only include the leading order term
    geff = 1.5
    omega = 2.**0.5*(1.-epsilon**2)**0.5
    width = 2.**0.5*epsilon

    period = twopi / omega
    if (epsilon > 1.e-5) then
       length = -log(rat/6./proamp)/width
    else
       length = 50.
    endif

! Check for some Light crossing criterion
!    if (length < 2.*period) length = 2.*period
    dx = length / nlat
    x0 = length / 2.
    dk = twopi / length

    print*,"period is ",period," length is ",length
  end subroutine set_coeffs

! Adjust this as needed, maybe use phi_max not v
  real(dl) function output_param()
    output_param = epsilon
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
!
  real(dl) function effective_mass(x,t)
    real(dl) :: x,t

    real(dl) :: p1, phitmp
    real(dl) :: xcur, tcur

    xcur = width*(x-x0)
    tcur = omega*t

    p1 = proamp / cosh(xcur)
!    phitmp = p1*cos(tcur) + geff*p1**2*(2.*cos(tcur) + cos(2.*tcur)-3.)
    phitmp = p1*cos(tcur) + geff*p1**2*(2.*cos(tcur) - 3.)
    effective_mass = 3.*(phitmp + 1.)**2 - 1.

  end function effective_mass

! Make sure that my choice of stencils for the |grad|^2 operator is actually the energy conserving choice (really only matters if I decide to include expansion)
  subroutine energy_density()

  end subroutine energy_density

end module Model
