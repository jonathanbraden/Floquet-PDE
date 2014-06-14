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
  real(dl), parameter :: alphamin = 20.
  real(dl) :: kval
  real(dl) :: xkink, vrun

  integer :: numsteps
  real(dl) :: tstep

  logical :: check_courant
  real(dl) :: tevolve1(2), tevolve2(2)

! Set parameters controlling integration
  outstep=256  ! total number of output per period
  numperiods=10

! Set model parameters
! To do, I need to convert a v value into an x0 for the kinks
  vrun = 0.01
  xkink = (1+vrun**2)**0.5*log(0.5*vrun)
  call set_coeffs(vrun)
!  kval = 4.006e-4  ! first stability band
!  kval = 0.065  ! large floquet exponent
!  kval = 4.925e-3 ! higher stability band
!  kval = 0.0006  ! 2nd instability band
!  kval = 0.0002
!  kval = 0.83  ! high order stability band
!  kval = 0.79   ! high order unstable mode
  kval = 0.05
  call set_k2eff(k2eff(kval))

! To do: In here put a loop to make sure we don't violate Courant
! All I'll need to do is keep doubling the number of time steps until the condition is met
  check_courant = .false.
  i=1
  print*, "dx = ",dx
  do while (.not.check_courant)
     numsteps = nlat*i
     tstep = period / dble(numsteps)
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
     delphi(i,1) = 1./cosh(xcur-xkink)
  enddo
  delphi(:,2) = 0.
  delphi(:,3) = 0.
  delphi(:,4) = delphi(:,1)

  do i=1,nlat
     write(98,'(2(F10.4,2x),6(ES14.7,2x))') (i*dx-length/2.), 0.,   &
          delphi(i,1), delphi(i,2), delphi(i,3), delphi(i,4), &
          effective_mass(i*dx,0.)
  enddo
  write(98,*) ""

  ! Now evolve my two orthogonal initial states
  call init_arrays()
  tevolve1 = 0.
  tevolve2 = 0.
  do j=1, outstep*numperiods
     tcur = dble(j*outstepsize)*tstep

     fld(1:nlat) = delphi(:,1)
     fldp(1:nlat) = delphi(:,2)
     call step(tstep, outstepsize, fld, fldp, tevolve1)
     delphi(:,1) = fld(1:nlat)
     delphi(:,2) = fldp(1:nlat)

     fld(1:nlat) = delphi(:,3)
     fldp(1:nlat) = delphi(:,4)
     call step(tstep, outstepsize, fld, fldp, tevolve2)  ! need to fix this
     delphi(:,3) = fld(1:nlat)
     delphi(:,4) = fldp(1:nlat)

     do i=1,nlat
        xcur = dble(i)*dx - length/2.
        write(98,'(2(F12.4,2x),6(ES14.7,2x))') xcur, tcur, &
             delphi(i,1), delphi(i,2), delphi(i,3), delphi(i,4), &
             effective_mass(i*dx,tcur)
     enddo
     write(98,*) ""
  enddo
end program evolve_initcond
