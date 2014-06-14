!
! Evolve the fluctuation equation from an initial state that is specified by hand
! We wish to use this to more deeply study the production of fluctuations during wall collisions, for which the initial state is taken to be the pair of orthogonal states corresponding to a bound state fluctuation.
!
! This is much simpler than the other driving programs in this directory
! since I never need to compute the fundamental matrix
!
program evolve_initcond
  use evolve

  implicit none

  integer :: i, j
  real(dl) :: delphi(1:nlat,4)  ! store the evolving fields
  integer :: outstep, numperiods, outstepsize
  real(dl) :: xcur, tcur
  real(dl) :: alphacheck
  real(dl), parameter :: alphamin = 5.
  real(dl) :: kval
  real(dl) :: rad_init

  integer :: numsteps
  real(dl) :: tstep

  logical :: check_courant
  real(dl) :: tevolve1(2), tevolve2(2)

! Set parameters controlling integration
  outstep=256  ! total number of output per period
  numperiods=5

! Set model parameters
! To do, I need to convert a v value into an x0 for the kinks
  rad_init = 10.
  call set_coeffs(rad_init)
  kval = 0.5   !0.51
  call set_k2eff(k2eff(kval))

! To do: In here put a loop to make sure we don't violate Courant
! All I'll need to do is keep doubling the number of time steps until the condition is met
  check_courant = .false.
  i=1
  print*, "dx = ",dx
  do while (.not.check_courant)
     numsteps = nlat*i
     tstep = get_period() / dble(numsteps)
     print*, tstep
     check_courant = (get_dx()/tstep .gt. alphamin)
     i=i*2
     print*,check_courant
  end do

  ! Now evolve the system forward in time
  if (outstep .gt. numsteps) numsteps = outstep
  outstepsize = numsteps / outstep
  if ( ( outstepsize*outstep) .ne. numsteps) then
     print*, "Error, outsteps don't divide total steps"
     stop
  endif

  open(unit=98, file="evolved_particles.dat")
  write(98,*) "# x  t   phi_1  pi_1   phi_2  pi_2  mass^2"

! Start with the field mode
  do i=1,nlat
     xcur = i*dx - length/2.
     delphi(i,1) = 2.**(-0.5)/cosh((xcur-rad_init)/2.**0.5)**2
  enddo
  delphi(:,2) = 0.
  delphi(:,3) = 0.
  delphi(:,4) = delphi(:,1)

  do i=1,nlat
     xcur = i*dx -length/2.
     write(98,'(2(F10.4,2x),6(ES14.7,2x))') xcur, 0.,   &
          delphi(i,1), delphi(i,2), delphi(i,3), delphi(i,4), effective_mass(i*dx,0.)
  enddo
  write(98,*) ""

  ! Now evolve my two orthogonal initial states
  call init_arrays()
  tevolve1 = 0.
  tevolve2 = 0.
  do j=1, outstep*numperiods
     tcur = dble(j*outstepsize)*tstep

!     print*, " before ", tevolve1
     fld(1:nlat) = delphi(:,1)
     fldp(1:nlat) = delphi(:,2)
     call step(tstep, outstepsize, fld, fldp, tevolve1)
     delphi(:,1) = fld(1:nlat)
     delphi(:,2) = fldp(1:nlat)

!     print*," and 2 is",tevolve2
     fld(1:nlat) = delphi(:,3)
     fldp(1:nlat) = delphi(:,4)
     call step(tstep, outstepsize, fld, fldp, tevolve2)  ! need to fix this
     delphi(:,3) = fld(1:nlat)
     delphi(:,4) = fldp(1:nlat)

!     print*, "tevolve 1 = ",tevolve1, " tevolve2 = ",tevolve2," tcur = ",tcur

     do i=1,nlat
        xcur = dble(i)*dx - length/2.
        write(98,'(2(F10.4,2x),6(ES14.7,2x))') xcur, tcur, &
             delphi(i,1), delphi(i,2), delphi(i,3), delphi(i,4), &
             effective_mass(i*dx,tcur)
     enddo
     write(98,*) ""
  enddo
end program evolve_initcond
