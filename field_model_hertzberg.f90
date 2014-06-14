module Model

implicit none

integer, parameter :: dl = kind(1.d0)
real(dl), parameter :: twopi = 6.283185307179586476925867665590_dl
!
! Model parameters, such as lattice size, etc
!
integer, parameter :: nlat = 512  ! Size of the lattice being used'
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

! Model specific (ie. sine-gordon) parameters
real(dl) :: xscl, tscl, amposc
real(dl) :: geff, epsilon
real(dl),parameter :: rat=1.e-6  ! value of background breather at boundaries

contains
!
! This is simply to speed up the computations.  We compute the values of the space-dependent coefficients once so we don't have to compute them many times
!
  subroutine set_coeffs(epsset, gset)
    real(dl) :: epsset, gset

    real(dl) :: width

    geff = gset
    epsilon = epsset

    amposc = 4.*epsilon*geff
    tscl = sqrt(1.-epsilon**2)
    xscl = epsilon

! amposc = 4.*epsilon*geff/(1-epsilon**2)
! xscl = epsilon/sqrt(1-epsilon**2)
! tscl = sqrt(1-epsilon)
! amp = geff
! k2 = 0.3**2/(1-epsilon**2)

    period = twopi/tscl  ! check this
! Set the length so that the field has damped down to some value R at the boundaries (needs to be set on a bg by bg basis
    length = 300.  !-2. * log(rat/amposc) / epsilon

    dx = length / nlat
    x0 = length / 2.
    dk = twopi / length

    print*,"period is ",period," length is ",length
  end subroutine set_coeffs

! Adjust this as needed, maybe use phi_max not v
  real(dl) function output_param()
    output_param = geff
  end function output_param

!
! To do, add in a subroutine to make a breather for actually running bg evolution
!

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

    real(dl) :: xeff, teff

    xeff = xscl * (x-x0)
    teff = tscl*t
    effective_mass = amposc*cos(teff)/cosh(xeff)
  end function effective_mass

! Make sure that my choice of stencils for the |grad|^2 operator is actually the energy conserving choice (really only matters if I decide to include expansion)
  subroutine energy_density()

  end subroutine energy_density

end module Model
