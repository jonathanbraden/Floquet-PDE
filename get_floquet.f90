program get_floquet

  use omp_lib
  use fundamental_matrix

  implicit none

  integer :: i, j, l
  logical :: check_courant
!  real*8, parameter :: delk = 0.005_dl
  real*8 :: delk
  real*8, parameter :: delv = 0.5
  real*8, parameter :: kmax = 12.5  !50.  ! maximum of k^2/m^2
  integer, parameter :: numk = 200   !200
  real*8, parameter :: deps = 1.
  real*8, parameter :: delamp = 0.01
  real*8, parameter :: alphamin = 20.
  real*8 :: kcur
  real*8 :: vtmp

  real*8, parameter :: pi = 3.141592653589793
  integer, parameter :: vlength = 3

  real*8, dimension(1:2*nlat,1:2*nlat) :: atmp

  call set_workarray(.false.)
  call init_basis()
  open(unit=99,file="floquet_chart.dat")
  write(99,*),'# nlat = ',nlat,' length = ',length,' tstep = ',numsteps
  write(99,*),'# amp   k^2    Floquet'

  do i=1,50
!
! Uncomment this for sine-Gordon
!
     ! uniformly space in 1/v
!     vtmp = 10.+i*1.
!     vtmp = 1./vtmp
!     vtmp = 0.5
!     vtmp = 1./(2.**0.5-1.)

!     call set_coeffs(vtmp)
     ! scan evenly over k^2T^2 for sine-gordon
!     delk = ( kmax / dble(numk) ) * vtmp**2 / (1.+vtmp**2)
!     call set_coeffs(5., 20.+i*0.5)
!     delk = 600. / dble(50) / get_period()**2

     call set_coeffs(i*0.01)
     delk = 1./numk

     check_courant = .false.
     l=1
     do while (.not.check_courant) 
        call set_numsteps(l*nlat)
        check_courant = (get_dx()/get_tstep() .gt. alphamin)
        l=l*2
     enddo

     do j=0,numk
        kcur = j*delk
        call set_k2eff(k2eff(kcur))
        call make_matrix()
        call get_exponents(99, kcur)
        print*,k2
        call flush(99)
     enddo
     write(99,*)
  enddo

!  do j=30,50
!     amp = dble(j)*delamp
!     amp = delamp*dble(j)
!     call set_coeffs(j*delvparam)
!     call set_coeffs(j*delamp, wid)
!     do i = -50,50
!        k2=dble(i)*dk
!        call make_matrix()
!        call get_exponents()
!     enddo
!     print*
!  enddo
end program get_floquet
