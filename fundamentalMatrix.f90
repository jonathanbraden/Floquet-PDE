module fundamentalMatrix
  use omp_lib
  use Model
  use evolve  ! fix this one up
  implicit none

! To do : encapsulate this information in a type
  real(dl), dimension(1:2*nlat,1:2*nlat) :: fmatrix
  real(dl), dimension(1:2*nlat) :: evalReal, evalImag
  real(dl), allocatable, dimension(:,:) :: leftVectors, rightVectors

! Storage for initial basis functions and inverse of initial matrix
  real(dl), dimension(1:nlat,1:nlat) :: basis
  real(dl), dimension(1:2*nlat,1:2*nlat) :: binv

! Declare work arrays, etc
  real(dl), dimension(:), allocatable :: workarray
  integer :: iworkarray

  type fundMatrix
     real(dl), allocatable, dimension(:,:) :: fmatrix
     real(dl), allocatable, dimension(:) :: evalReal, evalImag
     real(dl), allocatable, dimension(:,:) :: leftVector, rightVectors
  end type fundMatrix

  type workMatrix
     real(dl), allocatable, dimension(:) :: wArray
     integer :: iWArray
  end type workMatrix

contains

  subroutine createFundamentalMatrix(this, n)
    type(fundMatrix), intent(out) :: this
    allocate(fmatrix(n,n))
! Call subroutine to initialise work arrays here
  end subroutine createFundamentalMatrix

  subroutine makeFundamentalMatrix(model)
    type(fieldModel) :: model
    
    integer :: i, n
    type(myLattice) :: lat

    if (.not.set_step) print*, "Error, step size not set before evolving"
    print*,"Computing fundamental matrix"

!$OMP PARALLEL PRIVATE(lat,n)
    call makeLattice(lat)
    n = lat%nlat
!$OMP DO
    do i=1,nlat
       lat%fld(1:n) = basis(:,i)
       lat%fldp(1:n) = 0._dl
       call resetLatticeClock(lat)
       call step( ,numsteps)
       fmatrix(i,1:n) = fld(1:n)
       fmatrix(i,n+1:2*n) = fldp(1:n)

       lat%fld(:) = 0._dl
       lat%fldp(1:n) = basis(:,i)
       call resetLatticeClock(lat)
       call step( ,numsteps)
       fmatrix(n+i,1:n) = fld(1:n)
       fmatrix(n+i,n+1:2*n) = fldp(1:n)
    enddo
!$OMP END DO
!$OMP END PARALLEL

    fmatrix = matmul(binv,fmatrix)
  end subroutine makeFundamentalMatrix

  subroutine getExponents(this)
    type(fundMatrix) :: this
  end subroutine getExponents

  subroutine initialiseBasis(n)
    integer, intent(in) :: n

    integer :: i,j
    real(dl), dimension(2*n,2*n) :: atmp

! Compute inverse of initial basis matrix
    do i=1,n
       atmp(i,1:n) = basis(:,i)
       atmp(i,n+1:2*n) = 0._dl
       atmp(n+i,1:n) = 0._dl
       atmp(n+i,n+1:2*n) = basis(:,i)
    enddo
    call invert_matrix(atmp, binv, 2*n)
    
  end subroutine initializeBasis

  subroutine set_workarray(get_vectors)
    logical :: get_vectors
    integer :: ierror
    real(dl) :: dummy(1:1)

    if (.not.get_vectors) then
       call DGEEV('N','N', 2*nlat, fmatrix, 2*nlat, evalreal, evalimag, &
            dummy, 1, dummy, 1, dummy, -1, ierror)
    else
       allocate( leftvectors(1:2*nlat,1:2*nlat) )
       allocate( rightvectors(1:2*nlat,1:2*nlat) )
       call DGEEV('V','V', 2*nlat, fmatrix, 2*nlat, evalreal, evalimag, &
            leftvectors, 2*nlat, rightvectors, 2*nlat, dummy, -1, ierror)
    endif

    if (ierror == 0) then
       iworkarray = int(dummy(1))
       allocate(workarray(1:iworkarray))
    else
       print*, "Error allocating workspace array, exiting"
       stop
    endif
  end subroutine set_workarray

  

!
! To do, instead of killing program, just exit this subroutine if LU decomposition fails
! Put this into a separate module file
!
  subroutine invert_matrix(A,Ainv,n)
    real(dl), dimension(n,n), intent(in) :: A
    real(dl), dimension(n,n), intent(out) :: Ainv
    integer :: n
!!!!!!!!!!!!!!!!!!!!!!!!
! Workspace for LAPACK !
!!!!!!!!!!!!!!!!!!!!!!!!
    real(dl), allocatable, dimension(:) :: work
    integer :: lwork
    integer :: info, lda
    integer, allocatable, dimension(:) :: IPIV
    integer :: deallocatestatus
    external DGETRF
    external DGETRI

    lda = n
    lwork = n*n
    allocate(work(lwork))
    allocate(IPIV(n))

! Perform LU decomposition
    Ainv = A
    call DGETRF( n, n, Ainv, lda, IPIV, info )
    if (info.eq.0) then
       print*, "LU decomposition successful"
    elseif (info < 0) then
       print*, "LU decomposition: illegal value at ", inflo
       stop
    elseif (info > 0) then
       print*, "Singular Matrix U = 0 at ",info
    endif

    call DGETRI( n, Ainv, lda, IPIV, work, lwork, info )
    if (info.ne.0) then
       stop "Matrix inversion failder!"
    else
       print*, "Inverse Successful"
    endif
    
    deallocate(ipiv, stat=deallocatestatus)
    deallocate(work, stat=deallocatestatus)
  end subroutine invert_matrix

end module fundamentalMatrix
