module Model

implicit none

integer, parameter :: dl = kind(1.d0)
real(dl), parameter :: twopi = 6.283185307179586476925867665590_dl
!
! Model parameters, such as lattice size, etc
!
integer, parameter :: nlat = 1024 ! Size of the lattice being used'
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
real(dl) :: v
real(dl) :: gamma
real(dl),parameter :: rat=1.e-5  ! value of background breather at boundaries

contains
!
! This is simply to speed up the computations.  We compute the values of the space-dependent coefficients once so we don't have to compute them many times
!
  subroutine set_coeffs(vparam)
    real(dl) :: vparam

    v = vparam
    gamma = 1./(1+vparam**2)**0.5

    period = twopi / gamma / v
! Due to symmetries of the SG equation, only need half a period of the breather
! to obtain the full period of V''
    period = 0.5*period

! Set the length so that the field has damped down to some value R at the boundaries (needs to be set on a bg by bg basis
    length = -8.*2.*2.*log(rat*v/2.)/gamma 

! Check for light crossing criterion
!    if (length < period) length = period
    dx = length / nlat
    x0 = length / 2.
    dk = twopi / length

    print*,"period is ",period," length is ",length
  end subroutine set_coeffs

! Adjust this as needed, maybe use phi_max not v
  real(dl) function output_param()
    output_param = v
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

  subroutine init_background()
    call make_breather(fldbg(1:nlat), fldpbg(1:nlat))
  end subroutine init_background

  subroutine make_breather(fvals, fpvals)
    real(dl), dimension(1:nlat) :: fvals, fpvals
    real(dl) :: vval

    integer :: i
    real(dl) :: x, xtmp
    real(dl) :: ptemp

    do i=1,nlat
       x=i*dx
       ptemp = v*cosh(gamma*(x-x0))
       fvals(i) = 4.*atan(1/ptemp)
    enddo
    fpvals = 0.
  end subroutine make_breather

!
! Model specific definition of the effective mass
!
  real(dl) function effective_mass(x,t)
    real(dl) :: x,t
    real(dl) :: ptemp

    ptemp = v*cosh(gamma*(x-x0))
    effective_mass = cos(4. * atan(cos(gamma*v*t)/ptemp) )

  end function effective_mass

  real(dl) function field_bg(x,t)
    real(dl), intent(in) :: x,t
    real(dl) :: ptemp

    ptemp = v*cosh(gamma*(x-x0))
    field_bg = 4.*atan(cos(gamma*v*t)/ptemp)
  end function field_bg

  real(dl) function effective_mass_bg(phi)
    real(dl) :: phi
    effective_mass_bg = cos(phi)
  end function effective_mass_bg

  real(dl) function vprime(phi)
    real(dl) :: phi
    vprime = sin(phi)
  end function vprime

! Make sure that my choice of stencils for the |grad|^2 operator is actually the energy conserving choice (really only matters if I decide to include expansion)
  subroutine energy_density()

  end subroutine energy_density

end module Model
