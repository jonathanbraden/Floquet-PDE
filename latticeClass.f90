module myLattice
  implicit none
  integer, parameter :: dl = kind(1.d0)
  real(dl), parameter :: twopi = 6.283185307179586476925867665590_dl

  type lattice
     integer :: nlat, pad
     integer :: imin, imax
     real(dl), allocatable :: fld, fldp
     real(dl), dimension(1:2) :: tGlobal  ! change this to number of Hamiltonians
     real(dl) :: dt, dx, len, kmax
  end type lattice

contains

  subroutine createLattice(this,n,pad,dt,len)
    type(lattice), intent(out) :: this
    integer, intent(in) :: n, pad
    real(dl), intent(in) :: dt, len
    integer :: iMin, iMax

    iMin = 1-pad; iMax = n+pad
    this%nlat = n; this%pad = pad
    this%imin = iMin; this%imax = iMax
    allocate( this%fld(iMin:iMax),this%fldp(iMin:iMax) )
    tGlobal = 0.

    this%dt = dt
    this%len = len; this%dx = len / dble(n)
    this%kmax = twopi / len
  end subroutine createLattice

  subroutine destroyLattice(this)
    type(lattice), intent(inout) :: this

    deallocate(this%fld,this%fldp)
  end subroutine destroyLattice
end module myLattice
