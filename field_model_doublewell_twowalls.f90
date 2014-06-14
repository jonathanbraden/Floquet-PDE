module Model

implicit none

integer, parameter :: dl = kind(1.d0)
real(dl), parameter :: twopi = 6.283185307179586476925867665590_dl
!
! Model parameters, such as lattice size, etc
!
integer, parameter :: nlat = 64 ! Size of the lattice being used'
integer, parameter :: pad = 2  ! padding to implement boundary conditions
! implement this feature in a nicer way
integer, parameter :: imin=1-pad, imax=nlat+pad

! make these allocatable
real(dl), dimension(imin:imax) :: fld, fldp, fldbg, fldpbg
real(dl) :: rcur, rpcur
real(dl) :: period
real(dl) :: length
real(dl) :: dx
real(dl) :: x0
real(dl) :: dk

! Model specific (ie. sine-gordon) parameters
real(dl) :: rmax
real(dl), parameter :: width=2.**0.5
real(dl),parameter :: rat=1.e-6  ! value of background breather at boundaries

contains
!
! This is simply to speed up the computations.  We compute the values of the space-dependent coefficients once so we don't have to compute them many times
!
  subroutine set_coeffs(ri, per)
    real(dl), intent(in) :: ri, per

    rmax = ri
    period = per  ! compute this via a quadrature integration

    length = 2.*(rmax+8.*width)

! Check for light crossing criterion
!    if (length < period) length = period
    dx = length / nlat
    x0 = length / 2.
    dk = twopi / length

    print*,"period is ",period," length is ",length
  end subroutine set_coeffs

! Adjust this as needed, maybe use phi_max not v
  real(dl) function output_param()
    output_param = period
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

  real(dl) function effective_mass(x,t)
    real(dl), intent(in) :: x,t
    real(dl) :: xeff
    real(dl) :: rcur, rdotcur, gamma
    real(dl) :: phi_cur

    xeff = x-x0
! Analytic approximation for r
    rcur = 
    rdotcur = 
    gamma = (1.-rdotcur**2)**(-0.5)

    phi_cur = tanh() - tanh() - 1.
    effective_mass = 3.*phi_cur**2 - 1.
  end function effective_mass

!
! Model specific definition of the effective mass
! Temporary, use a fixed evolution for the background
!
  real(dl) function effective_mass_stupid(x,t)
    real(dl), intent(in) :: x,t
    real(dl) :: phitemp
    real(dl), parameter :: winv=1./width

    real(dl), parameter :: rmin=-0.5
    real(dl) :: tmod, xeff

    tmod = mod(t,period)
    xeff = x-x0
    rcur = rmax + (rmin - rmax)*exp(-(tmod/period-0.5)**2/0.05)

    phitemp = tanh(winv*(xeff+rcur)) - tanh(winv*(xeff-rcur)) - 1
    effective_mass = 3.*phitemp**2 - 1.
  end function effective_mass_stupid

! Functions for background evolution
  real(dl) function vprime(rcur)
    real(dl) :: rcur
    
    vprime = 0.
  end function vprime

  real(dl) function kin_factor(rcur)
    real(dl) :: rcur

    kin_factor = 0.
  end function kin_factor

! Make sure that my choice of stencils for the |grad|^2 operator is actually the energy conserving choice (really only matters if I decide to include expansion)
  subroutine energy_density()

  end subroutine energy_density

end module Model
