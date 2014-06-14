module Model

implicit none

integer, parameter :: dl = kind(1.d0)
!
! Model parameters, such as lattice size, etc
!
integer, parameter :: nlat = 128   ! Size of the lattice being used'
integer, parameter :: pad = 2  ! padding to implement boundary conditions
real(dl), parameter :: length = 32.
real(dl), parameter :: dx = length / dble(nlat)
real(dl), dimension(1-pad:nlat+pad) :: fld, fldp

real(dl) :: amp
real(dl) :: width

contains

  subroutine set_coeffs(ampset, widthset)
    real(dl) :: ampset, widthset

    amp = ampset
    width = widthset
  end subroutine set_coeffs

  subroutine boundary_conditions()
    integer :: j

! subroutine to implement periodic boundary conditions
    do j=1,pad
       fld(1-j) = fld(nlat + 1 - j)
       fld(nlat+j) = fld(j)
       fldp(1-j) = fldp(nlat + 1 - j)
       fldp(nlat+j) = fldp(j)
    enddo
  end subroutine boundary_conditions

  real(dl) function effective_mass(x,t)
    real(dl) :: x,t

    real(dl), parameter :: x0 = length/2.
    real(dl) :: ptemp

    ptemp = profile(x)

! Note: I changed the amplitude from 3 to 3/2 between making charts.  New ones are going to be labelled by omega2 in file name

    effective_mass = 1.5*(amp*ptemp*sin(t) + 1.)**2
  end function effective_mass

  real(dl) function profile(x)
    real(dl) :: x

    real(dl), parameter :: x0=length/2. !, width = 3.**2

    profile = exp(-(x-x0)**2/(2.*width**2))
  end function profile

! Make sure that my choice of stencils for the |grad|^2 operator is actually the energy conserving choice (really only matters if I decide to include expansion)
  subroutine energy_density()

  end subroutine energy_density

end module Model
