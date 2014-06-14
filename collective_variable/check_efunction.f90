!
! This is a program to check the eigenfunctions that are extracted from the efunction program.
! I comput the appropriate eigenfunction in there.
! I then initialize the field using that eigenfunction and proceed to evolve it for a single period (or multiple periods as a real check) using the symplectic integrator.
! At the end of the day, we check to make sure that it has indeed grown by the correct factor
!
program evolve_efunc

  use fundamental_matrix
  use Hamiltonian

  implicit none

  integer :: i, j
  integer, parameter :: numexp=2  ! number of floquet modes to get
  real(dl) :: curvals(1:4*numexp)
  real(dl) :: efunc2(1:nlat,4*numexp)
  real(dl) :: floq_exponent(numexp)
  integer :: outstep, numperiods, outstepsize
  real(dl) :: xcur, tcur
  real(dl) :: alphacheck
  real(dl), parameter :: alphamin = 5.
  real(dl) :: kval

  logical :: check_courant
  real(dl) :: tevolve(2)
  real(dl), dimension(imin:imax) :: fld, fldp
  real(dl) :: rcur

! Set parameters controlling integration
  outstep=256  ! total number of output per period
  numperiods=16

! Set model parameters
  call set_coeffs(2.5, -0.4)
  kval = 0.005
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

  print*,"time step is ", get_tstep()
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
  call get_eigenfunction(efunc2, floq_exponent, numexp)
  print*,"done exponent"

  ! Now that I've got my eigenfunctions to start with, evolve them in time
  if (outstep .gt. numsteps) numsteps = outstep
  outstepsize = numsteps / outstep
  if ( ( outstepsize*outstep) .ne. numsteps) then
     print*, "Error, outsteps don't divide total steps"
     stop
  endif

!!!!
!! To do : Add a loop here to evolve all of the extracted eigenfunctions
!!!!

  open(unit=98, file="evolved_efunc.dat")
  fld(1:nlat) = efunc2(:,1)
  fldp(1:nlat) = efunc2(:,2)
  laplace(1:nlat) = fld(1:nlat)
  call laplacian_1d(nlat, laplace, Fk, twopi/length, planf, planb)
  fld(0) = fld(nlat)
  fld(nlat+1) = fld(1)
  rcur = r_evolution(0.)
  do i=1,nlat
     write(98,'(2(F10.4,2x),6(ES14.7,2x))') (i*dx-length/2.), 0., rcur, efunc2(i,1), efunc2(i,2), &
          effective_mass_rad(i*dx,rcur), laplace(i), (fld(i+1)+fld(i-1)-2.*fld(i))/dx**2 !, &
!          (-laplace(i)*efunc2(i,1)**2+effective_mass(i*dx,0._dl)*efunc2(i,1))
  enddo
  write(98,*) ""


  ! Now that I have the eigenfunction, put it into the symplectic integrator and check
  tevolve = 0.
  rcur = rmax
  do j=1, outstep*numperiods
     tcur = dble(j*outstepsize)*tstep

     fld(1:nlat) = efunc2(:,1)
     fldp(1:nlat) = efunc2(:,2)
     call step(tstep, outstepsize, fld, fldp, tevolve)
     efunc2(:,1) = fld(1:nlat)
     efunc2(:,2) = fldp(1:nlat)

! Since it's assumed I'm doing a spectral evolution here, get the laplacian
     laplace = fld(1:nlat)
     call laplacian_1d(nlat, laplace, Fk, twopi/length, planf, planb)
! This part is only so I can check the finite difference approximation below
     fld(0) = fld(nlat)
     fld(nlat+1) = fld(1)
     rcur = r_evolution(tcur)
     do i=1,nlat
        xcur = dble(i)*dx - length/2.
        write(98,'(2(F8.4,2x),8(ES14.7,2x))') xcur, tcur, rcur, &
             efunc2(i,1), efunc2(i,2), effective_mass_rad(i*dx, rcur), laplace(i), (fld(i+1)+fld(i-1)-2.*fld(i))/dx**2 !,  &
!             (-laplace(i)*efunc2(i,1)**2+effective_mass(i*dx,0._dl)*efunc2(i,1))
     enddo
     write(98,*) ""
  enddo
end program evolve_efunc
