!
! Module to produce the fundamental matrix for my Floquet system with spatially dependent mass term
! I solve the differential equation many times (for different k^2 values), stick this into the fundamental matrix, then get eigenvalues, eigenvectors
!

! Test case: normal uniform background, since the solution there is known

!!!
! Choose initial basis for the fundamental matrix
!  Sine basis is used for Dirichlet I.C.s (in which case I should use a sine-transform or else Chebychev instead of Fourier for the spectral derivatives)
!!!
!#define DELTABASIS 1
#define FOURIERBASIS 1
!#define SINEBASIS 1

module fundamental_matrix

  use evolve
  use omp_lib

  implicit none

! Here we have to be careful to avoid violating the Courant condition !!
!  integer,parameter :: numsteps = 256*2
!  real(dl), parameter :: twopi = 6.2831853071795864769252867665590_dl
!  real(dl), parameter :: period = twopi/sqrt(1-0.05**2)
!  real(dl), parameter :: period = twopi
!  real(dl), parameter :: tstep = period/(dble(numsteps))

  integer :: numsteps
  real(dl) :: tstep
  logical :: set_step

  real(dl), dimension(1:2*nlat,1:2*nlat) :: fmatrix  ! The fundamental matrix
  real(dl), dimension(1:2*nlat) :: evalreal, evalimag
  real(dl), allocatable, dimension(:,:) :: leftvectors, rightvectors
!  real(dl), dimension(1:2*nlat,1:2*nlat) :: leftvectors, rightvectors

! Storage for initial basis functions and inverse of initial matrix
  real(dl), dimension(1:nlat,1:nlat) :: basis
  real(dl), dimension(1:2*nlat,1:2*nlat) :: binv

! Declare this here so I'm not reallocating it every time I run a new k-value
  real(dl), dimension(:), allocatable :: workarray
  integer :: iworkarray
  real(dl), dimension(2) :: tg

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

  subroutine make_matrix()
    integer :: i

    real(dl), dimension(1-pad:nlat+pad) :: fld, fldp

    if (.not.set_step) print*, "Error, step size not set before evolving"
    print*,"Computing fundamental matrix"

!    allocate(fld(1-pad:nlat+pad))
!    allocate(fldp(1-pad:nlat+pad))

! Ok, this breaks when I try spectral parallelization
! To fix I probably want to make a lattice class and make a copy of each of them
!$OMP PARALLEL DO PRIVATE(fld, fldp, tg, laplace, Fk)
    do i=1,nlat    ! loop over all initial basis vectors
       fld(1:nlat) = basis(:,i)
       fldp(:) = 0.
       ! Evolve given initial conditions over on oscillation period
       tg = 0.
       call step(tstep, numsteps, fld, fldp, tg)
       ! Store in the appropriate row of the matrix
       fmatrix(i,1:nlat) = fld(1:nlat)
       fmatrix(i,nlat+1:2*nlat) = fldp(1:nlat)
       ! Now evolve for the unit momentum solution
       fld(:) = 0.
       fldp(1:nlat) = basis(1:nlat,i)
       tg = 0.
       call step(tstep, numsteps, fld, fldp, tg)
       fmatrix(nlat+i,1:nlat) = fld(1:nlat)
       fmatrix(nlat+i,nlat+1:2*nlat) = fldp(1:nlat)
    enddo
!$OMP END PARALLEL DO

! Check ordering in here.  Replace with DGEMM call
    fmatrix = matmul(binv,fmatrix)
  end subroutine make_matrix

!
! Modify this subroutine to initialize the basis function of your choice
!
! For a periodic lattice with even number of lattice sites (and periodically extended x_i+N = x_i, sin(2piNx/2) is 0 at all of the lattice sites, hence is not a degree of freedom
!
  subroutine init_basis()
    integer :: i, j
    real(dl), dimension(2*nlat,2*nlat) :: atmp
! Fourier basis currently assumes nlat is even
#ifdef FOURIERBASIS
    real(dl) :: xcur, nfrac
    nfrac = twopi / nlat
    basis(:,1) = 1.
    do i=1,nlat/2-1
       do j=1,nlat
          xcur = (j-1)*nfrac
          basis(j,2*i) = sin(xcur*i)
          basis(j,2*i+1) = cos(xcur*i)
       enddo
    enddo
    do j=1,nlat
       xcur = (j-1)*nfrac
       basis(j,nlat) = cos(xcur*nlat/2)
    enddo
#endif
#ifdef SINEBASIS
    real(dl) :: scur, nfrac
    nfrac = twopi / 2. / nlat
    do i=1,nlat
       do j=1,nlat
          xcur = (j-1)*nfrac
          basis(j,i) = sin(xcur*i)
       enddo
    enddo
#endif

#ifdef DELTABASIS
    basis = 0.
    do i=1,nlat
       basis(i,i) = 1.
    enddo
#endif

! Compute inverse of initial matrix
    do i=1,nlat
       atmp(i,1:nlat) = basis(:,i)
       atmp(i,nlat+1:2*nlat) = 0.
       atmp(nlat+i,1:nlat) = 0.
       atmp(nlat+i,nlat+1:2*nlat) = basis(:,i)
    enddo
! Invert our initial fundamental matrix and store
    call invert_matrix(atmp, binv, 2*nlat)

    call init_arrays()

    print*,"initialization of fundamental matrix completed"
  end subroutine init_basis

!
! Subroutine to invert an NxN matrix
!
! A: real(dl) - nxn matrix to invert
! Ainv: real(dl) - nxn array to store the inverse of A
! n : size of arrays
!
  subroutine invert_matrix(A, Ainv, n)
    real(dl), dimension(n,n), intent(in) :: A
    real(dl), dimension(n,n), intent(out) :: Ainv
    integer :: n

    real(dl), allocatable, dimension(:) :: work
    integer :: lwork

    integer :: info, lda
    integer, allocatable, dimension(:) :: IPIV

    integer :: deallocatestatus
  
    external DGETRF  ! why are these external?
    external DGETRI

  ! check how much space LAPACK needs here
    lda = n  ! should really assign separately as input
    lwork = n*n
    allocate (work(lwork))
    allocate (IPIV(n))

! Perform LU decomposition
    Ainv = A
    call DGETRF( N, N, Ainv, LDA, IPIV, INFO )
! Add some calls to make sure the LU decomp worked
    if (info.eq.0) then
       print*, "LU decomposition successful"
    elseif (info.lt.0) then
       print*, "LU decomposition: illegal value at ", info
       stop
    elseif (info.gt.0) then
       print*, "Singular Matrix U = 0 at ",info
    endif

! USE LU decomposition to compute inverse
! To do: finish workspace query to determine size of work array
!  call DGETRI(N, Ainv, LDA, IPIV, WORK, -1, INFO)
!  lwork = work(1)
!  deallocate(work, 
!  allocate(work(lwork))
    call DGETRI(N, Ainv, LDA, IPIV, WORK, LWORK, INFO)

    if (info.ne.0) then
       stop "Matrix inversion failed!"
    else
       print*, "Inverse Successful"
    endif

    ! clean up temporary workspace
    deallocate(ipiv, STAT=deallocatestatus)
    deallocate(work, STAT=deallocatestatus)

  end subroutine invert_matrix

! Since LAPACK includes a workspace allocation call to determine the optimal size of the workspace array, just do this once before we start repeatedly diagonalizing the matrix
  subroutine set_workarray(get_vectors)
    logical :: get_vectors

    integer :: ierror
    real(dl) :: dummy(1:1)
    
    if (.not.get_vectors) then
       call DGEEV('N','N', 2*nlat, fmatrix, 2*nlat, evalreal, evalimag, dummy, 1, dummy, 1, dummy, -1, ierror)
    else
       allocate(leftvectors(1:2*nlat,1:2*nlat))
       allocate(rightvectors(1:2*nlat,1:2*nlat))
       call DGEEV('V','V', 2*nlat, fmatrix, 2*nlat, evalreal, &
            evalimag, leftvectors, 2*nlat, rightvectors, 2*nlat, dummy, -1, ierror)
    endif

    if (ierror .eq. 0) then
       iworkarray = int(dummy(1))
       allocate(workarray(1:iworkarray))
    else
       print*, "Error allocating workspace array, exiting"
       stop
    endif
  end subroutine set_workarray

  subroutine get_exponents(filenum, kout)
    integer :: filenum
    real(dl) :: kout

    integer :: maxind(1), i
    integer, parameter :: nfloquet=5
    real(dl) :: maxexp(1:nfloquet)
    integer,parameter :: asize = 2*nlat
    
    integer :: ierror
    real(dl) :: dummy(1,1)
    
    ! Call LAPACK routine to get eigenvalues
    ! Change first 2 N's to V's to get left/right eigenvectors
    call DGEEV('N','N', asize, fmatrix, asize, evalreal, evalimag, dummy, 1, dummy, 1, workarray, iworkarray, ierror)

! First, get the modulus of the eigenvalues (we don't really care about the phases for the moment
    evalreal = evalreal**2+evalimag**2
    do i=1,nfloquet
       maxexp(i) = maxval(evalreal)
       maxind = maxloc(evalreal)
       evalreal(maxind(1)) = 0.
    enddo
    write(filenum,'( 12(ES21.14,2x) )') output_param(), kout, log(maxexp(:))/2.

  end subroutine get_exponents

  subroutine get_eigenfunction(eigenfuncs, exponent)
    real(dl), dimension(1:nlat, 1:4), intent(out) :: eigenfuncs
    real(dl), intent(out) :: exponent

    integer :: i, maxind(1)
    real(dl) :: maxexp
    integer, parameter :: asize = 2*nlat
    real(dl) :: vectmp(1:nlat,4)

    integer :: ierror
    real(dl) :: dummy(1,1)
    
    call DGEEV('V','V', asize, fmatrix, asize, evalreal, evalimag, leftvectors, &
         asize, rightvectors, asize, workarray, iworkarray, ierror)

    evalreal = evalreal**2 + evalimag**2
    maxexp = maxval(evalreal)
    maxind = maxloc(evalreal)
    evalreal(maxind(1)) = 0.

    open(unit=99,file="eigenfunction.dat")
    write(99,*) "# F(x), pi(x), F(x), pi(x), Coeffs(x4) "

! Now extract the Fourier coefficients
    eigenfuncs(1:nlat,1) = leftvectors(1:nlat, maxind(1))
    vectmp(:,1) = matmul(eigenfuncs(:,1),basis)
    eigenfuncs(1:nlat,2) = leftvectors(nlat+1:2*nlat,maxind(1))
    vectmp(:,2) = matmul(eigenfuncs(:,2),basis)
    eigenfuncs(1:nlat,3) = rightvectors(1:nlat,maxind(1))
    vectmp(:,3) = matmul(eigenfuncs(:,3),basis)
    eigenfuncs(1:nlat,4) = rightvectors(nlat+1:2*nlat,maxind(1))
    vectmp(:,4) = matmul(eigenfuncs(:,4),basis)
    do i=1,nlat
       write(99,'(9(ES14.7,2x))') i*dx, eigenfuncs(i,1), eigenfuncs(i,2), &
            eigenfuncs(i,3), eigenfuncs(i,4), &
            vectmp(i,1), vectmp(i,2), &
            vectmp(i,3), vectmp(i,4)
    enddo

    print*, "# amp, k2, lyapanov exponent, imaginary piece"
    print*, output_param(), k2, log(maxexp)/2., evalimag(maxind(1))
    exponent = 0.5*log(maxexp)
  end subroutine get_eigenfunction

end module fundamental_matrix
