module Model

  implicit none
  integer, parameter :: dl = kind(1.d0) ! replace with a constants module
  real(dl), parameter :: twopi = 6.283185307179586476925867665590_dl

  type fieldModel
     integer :: nlat
     integer :: pad
     integer :: imin, imax
     real(dl) :: length, dx, dk

     real(dl), allocatable :: fldbg, fldpbg
     real(dl) :: period

     integer :: nParams
     real(dl), allocatable :: modParams
  end type fieldModel

contains

!
! Think about hardcoding the lattice size?
!
  subroutine createModel(this,n,pad,params)
    type(fieldModel), intent(out) :: this
    integer, intent(in) :: n,pad
    real(dl), dimension(:), intent(in) :: params
    integer :: imin, imax

    imin = 1 - pad; imax = nlat + pad

    this%nlat = n
    this%pad = pad
    this%imin = imin
    this%imax = imax

    allocate( this%fldbg(imin:imax), this%fldpbg(imin:imax) )

    this%nParams = size(params)
    allocate( this%modParams(1:this%nParams) )
    call setParams(this, params)
  end subroutine createModel

  subroutine destroyModel(this)
    type(fieldModel), intent(inout) :: this
    
    deallocate( this%modParams )
    deallocate( this%fldbg, this%fldpbg )
  end subroutine destroyModel

  subroutine setParams(this, params)
    type(fieldModel), intent(inout) :: this
  end subroutine setParams

  subroutine writeOutput(this)
    type(fieldModel), intent(in) :: this
  end subroutine writeOutput

  elemental real(dl) function effective_mass(this,x,t)
    type(fieldModel), intent(in) :: this
    real(dl), intent(in) :: x,t

    effective_mass = 0._dl
  end function effective_mass

end module Model
