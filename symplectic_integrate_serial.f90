module evolve

  use Hamiltonian

implicit none

contains

  subroutine step(dt, nsteps, fld, fldp, tglob)
    real(dl) :: dt
    integer :: nsteps
    real(dl), dimension(imin:imax) :: fld, fldp
    real(dl), dimension(1:2) :: tglob

    call symp4(dt, nsteps, fld, fldp, tglob)
  end subroutine step

! The basic integration consists of taking products of second order accurate time integrations.  Therefore, it's easiest to just include the second order time step as a separate subroutine
!
! To increase numerical efficiency, the factors for consecutive time steps are needed.  This is so that the last step of this product and the first step of the next one can be combined

  subroutine symp_o2step(dt, w1, w2)
    real(dl) :: dt, w1, w2

    integer :: i
!    integer :: n_Hamiltonian_pieces = 2  ! number of split terms in Hamiltonian

    do i=2,n_Hamiltonian_pieces-1
       call Hamiltonian_split(w1*dt/2._dl,i)
    enddo
    call Hamiltonian_split(w1*dt, n_Hamiltonian_pieces)
    do i=n_Hamiltonian_pieces-1,2,-1
       call Hamiltonian_Split(w1*dt/2._dl,i)
    enddo
    call Hamiltonian_Split((w1+w2)*dt/2._dl,1)
    return
  end subroutine symp_o2step

!
! Now define the various orders of integrators
! This is the entire meat of the method
!

  subroutine symp2(dt, nsteps)
    real(dl) :: dt
    integer :: nsteps

    integer :: j

    call Hamiltonian_Split(dt/2._dl,1)
    do j=1,nsteps-1
       call symp_o2step(dt,1._dl,1._dl)
    enddo
    call symp_o2step(dt,1._dl,0._dl)
  end subroutine symp2

!
! The fourth order integrator
!
  subroutine symp4(dt, nsteps)
    real(dl) :: dt
    integer :: nsteps

    integer :: j

    real(dl), parameter :: w1 = 1._dl/(2._dl-2._dl**(1._dl/3._dl))
    real(dl), parameter :: w0 = 1._dl - 2._dl*w1

    call Hamiltonian_Split(w1*dt/2._dl, 1)
    do j=1,nsteps
       call symp_o2step(dt,w1,w0)
       call symp_o2step(dt,w0,w1)
       if (j.eq.nsteps) then
          call symp_o2step(dt,w1,0._dl)
       else
          call symp_o2step(dt,w1,w1)
       endif
    enddo
  end subroutine symp4

  subroutine symp6(dt, nsteps)
    real(dl) :: dt
    integer :: nsteps

  end subroutine symp6

end module evolve
