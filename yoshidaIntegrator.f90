module yoshidaIntegrator

  use Hamiltonian
  use myLattice

implicit none

contains

  subroutine step(this, nsteps)
    type(lattice), intent(inout) :: this
    integer, intent(in) :: nsteps

    call symp6(this, nsteps)
  end subroutine step

! The basic integration consists of taking products of second order accurate time integrations.  Therefore, it's easiest to just include the second order time step as a separate subroutine
!
! To increase numerical efficiency, the factors for consecutive time steps are needed.  This is so that the last step of this product and the first step of the next one can be combined

  subroutine symp_o2step(this, dt, w1, w2)
    type(lattice), intent(inout) :: this
    real(dl), intent(in) :: dt, w1, w2

    integer :: i

    do i=2,n_Hamiltonian_pieces-1
       call Hamiltonian_split(w1*dt/2._dl,i, fld, fldp, tglob)
    enddo
    call Hamiltonian_split(w1*dt, n_Hamiltonian_pieces, fld, fldp, tglob )
    do i=n_Hamiltonian_pieces-1,2,-1
       call Hamiltonian_Split(w1*dt/2._dl,i, fld, fldp, tglob)
    enddo
    call Hamiltonian_Split((w1+w2)*dt/2._dl,1, fld, fldp, tglob)
    return
  end subroutine symp_o2step

!
! Now define the various orders of integrators
! This is the entire meat of the method
!

  subroutine symp2(this, nsteps)
    type(lattice), intent(inout) :: this
    integer, intent(in) :: nsteps
    
    real(dl) :: dt
    integer :: j

    dt = this%dt

    call Hamiltonian_Split(this,dt/2._dl,1)
    do j=1,nsteps-1
       call symp_o2step(this,dt,1._dl,1._dl)
    enddo
    call symp_o2step(this,dt,1._dl,0._dl)
  end subroutine symp2

!
! The fourth order integrator
!
  subroutine symp4(this, nsteps)
    type(lattice), intent(inout) :: this
    integer, intent(in) :: nsteps

    real(dl) :: dt
    integer :: j

    real(dl), parameter :: w1 = 1._dl/(2._dl-2._dl**(1._dl/3._dl))
    real(dl), parameter :: w0 = 1._dl - 2._dl*w1

    dt = this%dt

    call Hamiltonian_Split(w1*dt/2._dl, 1, fld, fldp, tglob)
    do j=1,nsteps
       call symp_o2step(dt,w1,w0, fld, fldp, tglob)
       call symp_o2step(dt,w0,w1, fld, fldp, tglob)
       if (j.eq.nsteps) then
          call symp_o2step(dt,w1,0._dl, fld, fldp, tglob)
       else
          call symp_o2step(dt,w1,w1, fld, fldp, tglob)
       endif
    enddo
  end subroutine symp4

  subroutine symp6(this, nsteps)
    type(lattice) :: this
    integer :: nsteps

    real(dl), parameter :: w1 = -1.17767998417887100694641568096431573_dl
    real(dl), parameter :: w2 = 0.235573213359358133684793182978534602_dl
    real(dl), parameter :: w3 = 0.784513610477557263819497633866349876
    real(dl), parameter :: w0 = 1._dl - 2._dl*(w1+w2+w3)

    real(dl) :: dt
    integer :: j

    dt = this%dt

    call Hamiltonian_Split(w3*dt/2._dl, 1)
    do j=1,nsteps
       call symp_o2step(this,dt,w3,w2)
       call symp_o2step(this,dt,w2,w1)
       call symp_o2step(this,dt,w1,w0)
       call symp_o2step(this,dt,w0,w1)
       call symp_o2step(this,dt,w1,w2)
       call symp_o2step(this, dt,w2,w3)
       if (j.eq.nsteps) then
          call symp_o2step(this,dt,w3,0._dl)
       else
          call symp_o2step(this,dt,w3,w3)
       endif
    enddo
  end subroutine symp6

end module yoshidaIntegrator
