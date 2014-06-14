!
! This is a program to check the eigenfunctions that are extracted from the efunction program.
! I comput the appropriate eigenfunction in there.
! I then initialize the field using that eigenfunction and proceed to evolve it for a single period (or multiple periods as a real check) using the symplectic integrator.
! At the end of the day, we check to make sure that it has indeed grown by the correct factor
!
program evolve_efunc

  use fundamental_matrix
  use evolve

  implicit none

  integer :: i, j
  real(dl) :: curvals(1:4)
  real(dl) :: efunc2(1:nlat,4)
  real(dl) :: floq_exponent
  integer :: outstep, numperiods, outstepsize
  real(dl) :: xcur, tcur
  real(dl) :: alphacheck
  real(dl), parameter :: alphamin = 20.
  real(dl) :: kval, vval

  logical :: check_courant
  real(dl) :: tevolve(2)

! Set parameters controlling integration
  outstep=128  ! total number of output per period
  numperiods=10

! Some samples for Sine-Gordon
  vval = 0.01
  kval = 0.05
!  vval = 1./(2.**0.5-1.)
!  kval = 0.2
!  vval = 1.
!  kval = 0.35
  call set_coeffs(vval)
!  kval = 52.*vval**2/(1.+vval**2)

! Samples for internal mode
!  call set_coeffs(0.5)
!  kval = 0.4

  call set_k2eff(k2eff(kval))
  
! To do: In here put a loop to make sure we don't violate Courant
! All I'll need to do is keep doubling the number of time steps until the condition is met
  check_courant = .false.
  i=1
  print*, "dx = ",dx
  do while (.not.check_courant)
     call set_numsteps(nlat*i)
     print*, tstep
     check_courant = (get_dx()/get_tstep() .gt. alphamin)
     i=i*2
     print*,check_courant
  end do

!  outstepsize = numsteps / outstep
!  if ( (outstepsize*outstep) .ne. numsteps) then
!     print*,"Error, outsteps don't divide total steps"
!     stop
!  endif

  call set_workarray(.true.)
  call init_basis()
  print*,"matrix"
  call make_matrix()
  print*,"matrix done"
  call get_eigenfunction(efunc2, floq_exponent)
  print*,"done exponent, floquet is ",floq_exponent

  ! Now that I've got my eigenfunctions to start with, evolve them in time
  if (outstep .gt. numsteps) numsteps = outstep
  outstepsize = numsteps / outstep
  if ( ( outstepsize*outstep) .ne. numsteps) then
     print*, "Error, outsteps don't divide total steps"
     stop
  endif

  open(unit=98, file="evolved_efunc.dat")
  fld(1:nlat) = efunc2(:,1)
  fldp(1:nlat) = efunc2(:,2)
  laplace(1:nlat) = fld(1:nlat)
  call laplacian_1d(nlat, laplace, Fk, twopi/length, planf, planb)
  fld(0) = fld(nlat)
  fld(nlat+1) = fld(1)
  do i=1,nlat
     write(98,'(2(ES21.14,2x),8(ES21.14,2x))') (i*dx-length/2.), 0., efunc2(i,1), efunc2(i,2), &
          effective_mass(i*dx,0._dl), laplace(i), (fld(i+1)+fld(i-1)-2.*fld(i))/dx**2, &
          (-laplace(i)*efunc2(i,1)**2+effective_mass(i*dx,0._dl)*efunc2(i,1))
  enddo
  write(98,*) ""


  ! Now that I have the eigenfunction, put it into the symplectic integrator and check
  tevolve = 0.
  do j=1, outstep*numperiods  !(numsteps*numperiods), outstepsize
     tcur = dble(j*outstepsize)*tstep
!     print*,tcur, j, tstep

     fld(1:nlat) = efunc2(:,1)
     fldp(1:nlat) = efunc2(:,2)
     call step(tstep, outstepsize, fld, fldp, tevolve)
     efunc2(:,1) = fld(1:nlat)
     efunc2(:,2) = fldp(1:nlat)
!     fld(1:nlat) = efunc2(:,3)
!     fldp(1:nlat) = efunc2(:,4)
!     call step(tstep, outstep)
!     efunc2(:,3) = fld(1:nlat)
!     efunc2(:,4) = fldp(1:nlat)

! Since it's assumed I'm doing a spectral evolution here, get the laplacian
     laplace = fld(1:nlat)
     call laplacian_1d(nlat, laplace, Fk, twopi/length, planf, planb)
! This part is only so I can check the finite difference approximation below
     fld(0) = fld(nlat)
     fld(nlat+1) = fld(1)

     do i=1,nlat
        xcur = dble(i)*dx - length/2.
        write(98,'(2(F12.4,2x),6(ES14.7,2x))') xcur, tcur, &
             efunc2(i,1), efunc2(i,2), effective_mass(i*dx, tcur), laplace(i), (fld(i+1)+fld(i-1)-2.*fld(i))/dx**2,  &
             (-laplace(i)*efunc2(i,1)**2+effective_mass(i*dx,0._dl)*efunc2(i,1))
     enddo
     write(98,*) ""
  enddo
end program evolve_efunc
