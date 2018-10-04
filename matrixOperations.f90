module matrixOperations
  integer, parameter :: 
  implicit none
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Invert a matrix.  This is done using an LU decomposition.
!
! Parameter :
!    A    - Input (square matrix to invert)
!    Ainv - Output (inverse of matrix A)
!    n    - Input (size of A.  A is an nxn matrix)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine invert_matrix(A,Ainv,n)
    real(dl), dimension(n,n), intent(in) :: A
    real(dl), dimension(n,n), intent(out) :: Ainv
    integer, intent(in) :: n

    real(dl), allocatable, dimension(:) :: work
    integer :: lwork
    integer :: info, lda
    integer, allocatable, dimension(:) :: ipiv
    integer :: deallocatestatus

    external DGETRF, DGETRI  ! why are these external

    lda = n
    lwork = n*n
    allocate(work(lwork))
    allocate(ipiv(n))

    Ainv = A
    call DGETRF( n, n, Ainv, lda, ipiv, info )
! Speed test the case statement versus if...elseif...elseif
    select case (info)
       case (0)
          print*, "LU decomposition successful"
       case (:-1)
          stop "LU decomposition: illegal value at ", info
       case(1:)
          stop "Singular Matrix U in LU decompsition at ", info
    end select
    
    call DGETRI( n, Ainv, lda, ipiv, work, lwork, info )

    if (info == 0) then
       print*,"Matrix inversion successful"
    else
       stop "Matrix inversion failed!"
    endif

    deallocate(ipiv, stat=deallocatestatus)
    deallocate(work, stat=deallocatestatus)
  end subroutine invert_matrix

  subroutine get_eigenvalues(A, eValues, n)
    real(dl), dimension(n,n), intent(in) :: A
    real(dl), dimension(n), intent(out) :: eValues
    integer, intent(in) :: n

    call DGEEV('N','N', )
    
  end subroutine get_eigenvalues

  subroutine get_eigenvectors(A,n)
    real(dl), dimension(n,n), intent(in) :: A
    
  end subroutine get_eigenvectors

end module matrixOperations
