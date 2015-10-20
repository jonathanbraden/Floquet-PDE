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
real(dl) :: amp
contains
!
! This is simply to speed up the computations.  We compute the values of the space-dependent coefficients once so we don't have to compute them many times
!
  subroutine set_coeffs(ampSet,spatialPeriods)
    real(dl), intent(in) :: ampSet
    integer, intent(in) :: spatialPeriods

    amp = ampSet
    length = twopi * spatialPeriods
    dx = length / dble(nlat)
    x0 = 0.5_dl * length
    dk = twopi / length

    period = twopi
    print*,"period is ",period," length is ",length
  end subroutine set_coeffs

! Adjust this as needed, maybe use phi_max not v
  real(dl) function output_param()
    output_param = amp
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

    real(dl) :: phitmp
    real(dl) :: xcur

    xcur = x-x0

    phitmp = amp * sin(xcur) * cos(t)
    effective_mass = phitmp

  end function effective_mass

end module Model
