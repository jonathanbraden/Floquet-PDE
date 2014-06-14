program get_floquet

  use fundamental_matrix

  implicit none

  integer :: i, j
!  real*8 :: dk = 0.02_dl
!  real*8 :: delamp = 0.1_dl
!  real*8 :: delvparam = 0.1_dl

!  real*8 :: delamp = 0.05_dl
!  real*8 :: wid

!  wid = 3.

  call set_workarray()
  print*,'# nlat = ',nlat,' length = ',length,' tstep = ',numsteps
  print*,'# amp   k^2    Floquet'

  k2=0.

  do i=0,15
     call set_coeffs(0.05_dl,i*0.05_dl)
     k2 = 0.3**2/(1-0.05**2)
     call make_matrix()
     call get_exponents()
     print*
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
