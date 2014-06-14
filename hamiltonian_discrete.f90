module Hamiltonian

  use Model

  implicit none

  integer, parameter :: n_Hamiltonian_pieces = 2  ! number of terms in split Hamiltonian

  real(dl), dimension(1:n_Hamiltonian_pieces) :: tglobal   ! stores current value of the time that the fields and momenta are being stored at
  real(dl) :: k2

 contains

   subroutine set_k2eff(keff)
     real(dl) :: keff
     k2 = keff
   end subroutine set_k2eff

   subroutine Hamiltonian_Split(dt, term_index, fld, fldp, tglob)
     real(dl) :: dt
     integer :: term_index  ! labels which term of the split hamiltonian to use
     real(dl), dimension(imin:imax) :: fld, fldp
     real(dl), dimension(1:2) :: tglob

     select case(term_index)
     case(1)
        tglob(1) = tglob(1) + dt
        call Hamiltonian_field(dt, fld, fldp)
     case(2)
        tglob(2) = tglob(2) + dt
        call Hamiltonian_momentum(dt, fld, fldp, tglob(1))
     case default
        print*, "Undefined Hamiltonian term"
        stop
     end select
   end subroutine Hamiltonian_Split

   subroutine Hamiltonian_field(dt, fld, fldp)
     real(dl) :: dt
     real(dl), dimension(imin:imax) :: fld, fldp
     integer :: i

     do i=1,nlat
        fld(i) = fld(i) + fldp(i)*dt
     enddo
   end subroutine Hamiltonian_field

! 
! Just need to fill in this subroutine
!
! Should use a macro to define my laplacian, and place the potential somewhere
! Also, i want to scan of k^2 values, so this needs to be permitted somehow
!
   subroutine Hamiltonian_momentum(dt, fld, fldp, time)
     real(dl) :: dt
     real(dl), dimension(imin:imax) :: fld, fldp
     real(dl) :: time

     real(dl) :: lap
     integer :: i

     call boundary_conditions()

     do i=1,nlat
        lap = - (fld(i+2) + fld(i-2)) + 16.*(fld(i+1) + fld(i-1)) - 30.*fld(i)
        lap = lap / (12.*dx**2)

        fldp(i) = fldp(i) - dt * ( (k2 + effective_mass(i*dx,time) )*fld(i) - lap )
     enddo
   end subroutine Hamiltonian_momentum

end module Hamiltonian
