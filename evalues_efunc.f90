!
! Module to produce the fundamental matrix for my Floquet system with spatially dependent mass term
! I solve the differential equation many times (for different k^2 values), stick this into the fundamental matrix, then get eigenvalues, eigenvectors
!

! Test case: normal uniform background, since the solution there is known

module eigenfunctions

  use evolve

! Here we have to be careful to avoid violating the Courant condition !!
!  integer,parameter :: numsteps = 256*32 !*2**6 !256
  real(dl) :: tstep
  logical :: set_step
  integer :: numsteps

  ! use Lapack  ! find a good matrix diagonalization package to use

  real(dl), dimension(1:2*nlat,1:2*nlat) :: fmatrix  ! The fundamental matrix
  real(dl), dimension(1:2*nlat) :: evalreal, evalimag
  real(dl), dimension(1:2*nlat,1:2*nlat) :: leftvectors, rightvectors

! Declare this here so I'm not reallocating it every time I run a new k-value
! This depends on the number of lattice sites, so it needs to be adjusted if that changes
  real(dl), dimension(:), allocatable :: workarray
  integer :: iworkarray

! Add in the storage (matrix?) for the eigenvectors if I want them
! Alternatively, I could just find the eigenvector for the max eigenvalue of something

! To use LAPACK, I'll need some temporary work arrays

contains

  subroutine set_numsteps(nsteps)
    integer :: nsteps

    numsteps = nsteps
    tstep = period / dble(numsteps)
    set_step = .true.
  end subroutine set_numsteps

  integer function get_numsteps()
    get_numsteps = numsteps
  end function get_numsteps

  real(dl) function get_tstep()
    get_tstep = tstep
  end function get_tstep

!
! This includes the equivalent of a numerical double check to make sure everything is working.
! If I choose the initial time fundamental matrix to be something other than the identity matrix, then I must multiply by the inverse of the initial fundamental matrix on the left before I go ahead and calculate all the eigenvalues and eigenvectors (otherwise the result will be incorrect).  C.f. lectures on Floquet theory.  I've only done the most trivial test with multiplies by -1's on the identity
!
  subroutine make_matrix
    integer :: i

    if (.not.set_step) print*, "Error, step size not set before evolving"

    do i=1,nlat
       fld(:) = 0.
       fldp(:) = 0.
       fld(i) = 1.
       ! Evolve given initial conditions over on oscillation period
       tglobal = 0.
       call step(tstep, numsteps)
       ! Store in the appropriate row of the matrix
       fmatrix(i,1:nlat) = fld(1:nlat)
       fmatrix(i,nlat+1:2*nlat) = fldp(1:nlat)
       ! Now evolve for the unit momentum solution
       fld(:) = 0.
       fldp(:) = 0.
       fldp(i) = 1.
       tglobal = 0.
       call step(tstep, numsteps)
       fmatrix(nlat+i,1:nlat) = fld(1:nlat)
       fmatrix(nlat+i,nlat+1:2*nlat) = fldp(1:nlat)
    enddo

    set_step = .false.
  end subroutine make_matrix

! Since LAPACK includes a workspace allocation call to determine the optimal size of the workspace array, just do this once before we start repeatedly diagonalizing the matrix
  subroutine set_workarray
    integer :: ierror
    real(dl) :: dummy(1:1)

    call DGEEV('V','V', 2*nlat, fmatrix, 2*nlat, evalreal, evalimag, leftvectors, 2*nlat, rightvectors, 2*nlat, dummy, -1, ierror)
    ! finish filling this in when I start computing eigenvectors
!    call DGEEV('V','V', 2*nlat, fmatrix, 2*nlat, evalreal, evalimag, 

    if (ierror .eq. 0) then
       iworkarray = int(dummy(1))
       allocate(workarray(1:iworkarray))
    else
       print*, "Error allocating workspace array, exiting"
       stop
    endif
  end subroutine set_workarray

!
! To do, allow for obtaining more than just the largest eigenvalue
!
  subroutine get_exponents(eigenfuncs)
    real(dl), dimension(1:nlat,1:4), intent(out) :: eigenfuncs

    integer :: i, j, maxind(1)
    integer, parameter :: nfloquet=1
!    integer, dimension(1:nfloquet) :: maxind
    real(dl) :: maxexp(1:nfloquet)
    integer,parameter :: asize = 2*nlat
    
    integer :: ierror
    real(dl) :: dummy(1,1)

    ! All I need to do in this subroutine is get the eigenvalues
    ! In the future I should also get eigenvalues (and I should really consider the log of the matrix)
    
    ! Call LAPACK routine to get eigenvalues
    ! Change first 2 N's to V's to get left/right eigenvectors
    call DGEEV('V','V', asize, fmatrix, asize, evalreal, evalimag, leftvectors, &
         asize, rightvectors, asize, workarray, iworkarray, ierror)

    !call DGEEV('V','V', asize, fmatrix, asize, evalreal, evalimag, leftvectors, asize,rightvectors, asize, workarray, iworkarray, ierror)
! to do, call a workspace query (ie. use iworkarray = -1) to get the size of the work array

! First, get the modulus of the eigenvalues (we don't really care about the phases for the moment
    evalreal = evalreal**2+evalimag**2
    do i=1,nfloquet
       maxexp(i) = maxval(evalreal)
       maxind = maxloc(evalreal)
       evalreal(maxind(1)) = 0.

       open(unit=99,file="eigenfunction.dat")

       if (i.eq.1) then
          eigenfuncs(1:nlat,1) = leftvectors(1:nlat,maxind(1))
          eigenfuncs(1:nlat,2) = leftvectors(nlat+1:2*nlat, maxind(1))
          eigenfuncs(1:nlat,3) = rightvectors(1:nlat, maxind(1))
          eigenfuncs(1:nlat,4) = rightvectors(nlat+1:2*nlat, maxind(1))
          do j=1,nlat
!             write(99,'(F8.4,4(ES14.7,2x))') eigenfuncs(i,1), eigenfuncs(i,2)
             write(99,'(4(ES14.7,2x))') leftvectors(j,maxind), leftvectors(nlat+j,maxind),  &
                  rightvectors(j,maxind), rightvectors(nlat+j,maxind)
          enddo
       endif
    enddo
    close(99)
    print*,'# amp, k2, lyapanov exponent, imaginary piece'
    print*, amp, k2, log(maxexp(:))/(2.*period), evalimag(maxind(1))



! Compute the maximual Floquet exponent
!   do i=2,2*nlat
!      curexp = real(log(eigenvalue(i)))
!      if (curexp .gt. maxexp) then
!         maxind=i
!         maxexp = curexp
!      endif
!   enddo

! Here, if I want to get the eigenstate associated with the maximal eigenvalue, I should get it

  end subroutine get_exponents

end module eigenfunctions
