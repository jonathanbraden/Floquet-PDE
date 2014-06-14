module Model

implicit none

integer, parameter :: dl = kind(1.d0)
real(dl), parameter :: twopi = 6.283185307179586476925867665590_dl
!
! Model parameters, such as lattice size, etc
!
integer, parameter :: nlat = 128   ! Size of the lattice being used'
integer, parameter :: pad = 2  ! padding to implement boundary conditions

!real(dl), parameter :: length = 256.
!real(dl), parameter :: dx = length / dble(nlat)
real(dl), dimension(1-pad:nlat+pad) :: fld, fldp, fldbg, fldpbg

real(dl) :: period
real(dl) :: length
real(dl) :: dx
real(dl) :: x0  ! midpoint of the lattice

! Model specific parameters
real(dl) :: epsilon
real(dl), parameter :: rat=1.e-6  ! value of bg at boundaries
real(dl) :: proamp
real(dl) :: omega
real(dl) :: width  ! depends on choice of coords
! Next order in expansion
real(dl) :: g_2

contains

  subroutine set_coeffs(epsset)
    real(dl) :: epsset

    epsilon = epsset
    proamp = 2.*epsilon / sqrt(3.)
    g_2 = 1.5  ! check this
    omega = 2.**0.5*(1.-epsilon**2)**0.5
    width = 2.**0.5*epsilon

    length = -2.*log(rat/2.)/width  ! fix this to depend on epsilon!!
    period = twopi/omega   ! adjust based on coordinate system
    dx = length / dble(nlat)
    x0 = length / 2.

  end subroutine set_coeffs

  real(dl) function output_param()
    output_param = epsilon
  end function output_param

  real(dl) function get_dx()
    get_dx = dx
  end function get_dx

  real(dl) function get_period()
    get_period = period
  end function get_period

  real(dl) function k2eff(kreal)
    real(dl) :: kreal
    k2eff = kreal
  end function k2eff

!
! Change this subroutine to be make background or something
! Or add what's basically a dummy routine that calls this
! This will allow the driver program to not have to do anything
! other than call the same make_background subroutine regardless of the model
!  
  subroutine make_oscillon(fvals, fpvals)
    real(dl), dimension(1:nlat) :: fvals, fpvals

    integer :: i
    real(dl) :: x, xtmp

    do i=1,nlat
       x = i*dx
       fvals(i) = profile(x, 0._dl)
    enddo

    fvals = fvals + 1.
    fpvals = 0.

  end subroutine make_oscillon

  subroutine boundary_conditions()
    integer :: j

! subroutine to implement periodic boundary conditions
    do j=1,pad
       fld(1-j) = fld(nlat + 1 - j)
       fld(nlat+j) = fld(j)
       fldp(1-j) = fldp(nlat + 1 - j)
       fldp(nlat+j) = fldp(j)
    enddo

    do j=1,pad
       fldbg(1-j) = fldbg(nlat+1-j)
       fldbg(nlat+j) = fldbg(j)
    enddo

  end subroutine boundary_conditions

  real(dl) function effective_mass(x,t)
    real(dl) :: x,t

    real(dl) :: phitmp

    phitmp = profile(x,t)

    effective_mass = 3.*(phitmp + 1.)**2 - 1.
  end function effective_mass

  real(dl) function effective_mass_bg(phi)
    real(dl) :: phi
    effective_mass_bg = 3.*phi**2 - 1.
  end function effective_mass_bg

  real(dl) function vprime(phi)
    real(dl) :: phi
    vprime = phi*(phi**2 - 1._dl)
  end function vprime

  real(dl) function profile(xcur, tcur)
    real(dl) :: xcur
    real(dl) :: tcur
    real(dl) :: xcentered
    real(dl) :: p1

    xcentered = xcur - x0
    p1 = proamp / cosh(width*xcentered)
! First line is leading order in epsilon, 2nd line is next order    
    profile = p1*cos(omega*tcur) + g_2*p1**2*(cos(2.*omega*tcur)-3.) / 6.
  end function profile

! Make sure that my choice of stencils for the |grad|^2 operator is actually the energy conserving choice (really only matters if I decide to include expansion)
  subroutine energy_density()

  end subroutine energy_density

end module Model
